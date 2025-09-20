/*
 * Copyright Lamont Granquist, Sebastien Gaggini and the MechJeb contributors
 * SPDX-License-Identifier: LicenseRef-PD-hp OR Unlicense OR CC0-1.0 OR 0BSD OR MIT-0 OR MIT OR LGPL-2.1+
 */

using System;
using System.Collections.Generic;
using MechJebLib.Functions;
using MechJebLib.Primitives;
using MechJebLib.Utils;

namespace MechJebLib.PDG
{
    public partial class Ascent
    {
        private readonly AscentBuilder _input;
        private          double        _eccT;
        private          double        _gammaT;
        private          Optimizer?    _optimizer;
        private          double        _smaT;

        private double _vT;

        private Ascent(AscentBuilder builder)
        {
            _input = builder;
        }

        private List<Phase> _phases              => _input._phases;
        private V3          _r0                  => _input._r0;
        private V3          _v0                  => _input._v0;
        private V3          _u0                  => _input._u0;
        private double      _t0                  => _input._t0;
        private double      _mu                  => _input._mu;
        private double      _rbody               => _input._rbody;
        private double      _peR                 => _input._peR;
        private double      _apR                 => _input._apR;
        private double      _attR                => _input._attR;
        private double      _incT                => _input._incT;
        private double      _lanT                => _input._lanT;
        private double      _fpaT                => _input._fpaT;
        private double      _hT                  => _input._hT;
        private bool        _attachAltFlag       => _input._attachAltFlag;
        private bool        _lanflag             => _input._lanflag;
        private bool        _fixedBurnTime       => _input._fixedBurnTime;
        private Solution?   _solution            => _input._solution;
        private int         _optimizedPhase      => _input._optimizedPhase;
        private int         _optimizedCoastPhase => _input._optimizedCoastPhase;

        public void Run()
        {
            (_smaT, _eccT) = Astro.SmaEccFromApsides(_peR, _apR);

            using Optimizer.OptimizerBuilder builder = Optimizer.Builder()
                .Initial(_r0, _v0, _u0, _t0, _mu, _rbody)
                .TerminalConditions(_hT);

            if (_solution == null)
            {
                _optimizer = _fixedBurnTime
                    ? InitialBootstrappingFixed(builder)
                    : InitialBootstrappingOptimized(builder);
            }
            else
            {
                _optimizer = ConvergedOptimization(builder, _solution);
            }
        }

        public Optimizer? GetOptimizer() => _optimizer;

        private Optimizer ConvergedOptimization(Optimizer.OptimizerBuilder builder, Solution solution)
        {
            if (_fixedBurnTime)
            {
                ApplyEnergy(builder);
            }
            else
            {
                if (_attachAltFlag || _eccT < 1e-4)
                    ApplyFPA(builder);
                else
                    ApplyKepler(builder);
            }

            ApplyOldBurnTimesToPhases(solution);
            Optimizer pdg = builder.Build(_phases);
            using Solution? solution2 = pdg.Run(solution);

            if (!pdg.Success() || solution2 == null)
                throw new Exception("converged optimizer failed");

            (V3 rf, V3 vf) = solution2.TerminalStateVectors();

            (_, _, _, _, _, double tanof, _) =
                Astro.KeplerianFromStateVectors(_mu, rf, vf);

            if (_attachAltFlag || _fixedBurnTime || Math.Abs(Statics.ClampPi(tanof)) < Math.PI / 2.0)
                return pdg;

            ApplyFPA(builder);
            ApplyOldBurnTimesToPhases(solution2);

            using Optimizer pdg2 = builder.Build(_phases);
            using Solution? solution3 = pdg2.Run(solution2);

            return pdg2.Success() ? pdg2 : pdg;
        }

        private void ApplyFPA(Optimizer.OptimizerBuilder builder)
        {
            // If _attachAltFlag is NOT set then we are bootstrapping with ApplyFPA prior to
            // trying free attachment with Kepler and attR is invalid and we need to fix to
            // the PeR.  This should be fixed in the AscentBuilder by having more APIs than
            // just "SetTarget" that fixes this correctly there.
            double attR = _attachAltFlag ? _attR : _peR;

            (_vT, _gammaT) = Astro.ConvertApsidesTargetToFPA(_peR, _apR, attR, _mu);
            if (_lanflag)
                builder.TerminalFPA5(attR, _vT, _gammaT, _incT, _lanT);
            else
                builder.TerminalFPA4(attR, _vT, _gammaT, _incT);
        }

        private void ApplyKepler(Optimizer.OptimizerBuilder builder)
        {
            (_smaT, _eccT) = Astro.SmaEccFromApsides(_peR, _apR);

            if (_lanflag)
                builder.TerminalKepler4(_smaT, _eccT, _incT, _lanT);
            else
                builder.TerminalKepler3(_smaT, _eccT, _incT);
        }

        private void ApplyEnergy(Optimizer.OptimizerBuilder builder)
        {
            if (_lanflag)
                builder.TerminalEnergy4(_attR, _incT, _lanT);
            else
                builder.TerminalEnergy3(_attR, _fpaT, _incT);
        }

        private Optimizer InitialBootstrappingFixed(Optimizer.OptimizerBuilder builder)
        {
            ApplyEnergy(builder);

            List<Phase> bootphases = DupPhases(_phases);

            // FIXME: we may want to convert this to an optimized burntime circular orbit problem with an infinite upper stage for bootstrapping
            foreach (Phase p in bootphases)
                p.Unguided = false;

            using Optimizer pdg       = builder.Build(bootphases);
            pdg.N = 6;
            using Solution  solution  = pdg.InitialGuess(_incT);
            using Solution? solution2 = pdg.Run(solution);

            if (!pdg.Success() || solution2 == null)
                throw new Exception("Target unreachable (fixed bootstrapping)");

            ApplyOldBurnTimesToPhases(solution2);

            List<Phase> bootphases2 = DupPhases(_phases);

            Optimizer pdg2 = builder.Build(bootphases2);
            using Solution? solution3 = pdg2.Run(solution2);

            if (!pdg2.Success() || solution3 == null)
                throw new Exception("Target unreachable");

            return pdg2;
        }

        private Optimizer InitialBootstrappingOptimized(Optimizer.OptimizerBuilder builder)
        {
            /*
             * Initial bootstrapping with infinite upper stage and no coast, forced FPA attachment
             */

            ApplyFPA(builder);
            List<Phase> bootphases = DupPhases(_phases);

            // switch the optimized phase to the top stage of the rocket
            bootphases[_optimizedPhase].OptimizeTime = false;

            // set the top stage to infinite + optimized + guided
            bootphases[bootphases.Count - 1].Infinite     = true;
            bootphases[bootphases.Count - 1].Unguided     = false;
            bootphases[bootphases.Count - 1].OptimizeTime = true;

            Statics.DebugPrint("*** PHASE 1: DOING INITIAL INFINITE UPPER STAGE ***");
            using Optimizer pdg           = builder.Build(bootphases);
            pdg.N = 6;
            using Solution  solutionGuess = pdg.InitialGuess(_incT);
            using Solution? solution      = pdg.Run(solutionGuess);

            if (!pdg.Success() || solution == null)
                throw new Exception("Target unreachable (infinite ISP)");

            ApplyOldBurnTimesToPhases(solution);

            Statics.DebugPrint("*** PHASE 2: RESETTING ROCKET AND DOING PERIAPSIS ATTACHMENT ***");
            Optimizer pdg2 = builder.Build(_phases);
            using Solution? solution2 = pdg2.Run(solution);

            if (!pdg2.Success() || solution2 == null)
                throw new Exception("Target unreachable");

            if (_attachAltFlag || _eccT < 1e-4 )
                return pdg2;

            /*
             * relaxing to free attachment
             */

            ApplyKepler(builder);
            ApplyOldBurnTimesToPhases(solution2);

            Statics.DebugPrint("*** PHASE 3: RELAXING TO FREE ATTACHMENT ***");
            Optimizer pdg3 = builder.Build(_phases);
            using Solution? solution3 = pdg3.Run(solution2);

            if (!pdg3.Success() || solution3 == null)
            {
                Statics.DebugPrint("*** FREE ATTACHMENT FAILED, FALLING BACK TO PERIAPSIS ***");
                return pdg2;
            }

            // this should catch if free attachment picked the apoapsis accidentally
            if (solution2.Vgo(solution3.T0) < solution3.Vgo(solution3.T0))
            {
                Statics.DebugPrint($"*** PERIAPSIS ATTACHMENT IS MORE OPTIMAL ({solution2.Vgo(solution3.T0)} < {solution3.Vgo(solution3.T0)}) THAN FREE ATTACHMENT SOLN ***");
                return pdg2;
            }

            // this catches issues with circular orbits where the kepler constraints fail to converge well
            if (pdg2.PrimalFeasibility * 100 < pdg3.PrimalFeasibility)
            {
                Statics.DebugPrint( $"*** FREE ATTACHMENT PRIMAL FEASIBILITY ({pdg3.PrimalFeasibility}) IS MUCH WORSE THAN PERIAPSIS ({pdg2.PrimalFeasibility}) ***");
                return pdg2;
            }

            return pdg3;
        }

        private void ApplyOldBurnTimesToPhases(Solution oldSolution)
        {
            for (int i = 0; i < _phases.Count; i++)
            {
                if (!_phases[i].OptimizeTime)
                    continue;

                for (int j = 0; j < oldSolution.Segments; j++)
                {
                    if (!oldSolution.OptimizeTime(j))
                        continue;

                    if (oldSolution.CoastPhase(j) != _phases[i].Coast)
                        continue;

                    _phases[i].bt = Math.Min(oldSolution.Bt(j, _t0), _phases[i].tau * 0.99);
                }
            }
        }

        // FIXME: this obviously creates garbage and needs to return an IDisposable wrapping List<Phases> that has a pool
        private List<Phase> DupPhases(List<Phase> oldphases)
        {
            var newphases = new List<Phase>();

            foreach (Phase phase in oldphases)
                newphases.Add(phase.DeepCopy());

            return newphases;
        }

        public static AscentBuilder Builder() => new AscentBuilder();
    }
}
