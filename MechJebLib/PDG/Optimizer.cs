/*
 * Copyright Lamont Granquist, Sebastien Gaggini and the MechJeb contributors
 * SPDX-License-Identifier: LicenseRef-PD-hp OR Unlicense OR CC0-1.0 OR 0BSD OR MIT-0 OR MIT OR LGPL-2.1+
 */

using System;
using System.Collections.Generic;
using System.Threading;
using MechJebLib.Functions;
using MechJebLib.ODE;
using MechJebLib.Primitives;
using MechJebLib.Utils;
using static MechJebLib.Utils.Statics;
using static System.Math;

// TODO: unguided upper stages.
// TODO: inequality constraints for coast phase "burn" time.
namespace MechJebLib.PDG
{
    public partial class Optimizer : IDisposable
    {
        public enum OptimStatus { CREATED, SUCCESS, CANCELLED, FAILED }

        // XXX: why do phases hang off the optimizer and not off the problem?
        private readonly List<Phase> _phases;
        private          int         _numPhases => _phases.Count;

        private readonly Problem                  _problem;
        private readonly alglib.minnlcreport      _rep = new alglib.minnlcreport();
        private readonly alglib.ndimensional_sjac _constraintHandle;
        private          alglib.minnlcstate       _state = new alglib.minnlcstate();

        private CancellationToken _timeoutToken;
        public  int               Iterations;
        public  OptimStatus       Status;
        public  int               TerminationType;
        public  double            PrimalFeasibility;
        public  double            Cost;

        public Solution? Solution;

        private int K => 2 * N - 1;

        private Optimizer(Problem problem, IEnumerable<Phase> phases)
        {
            _phases           = new List<Phase>(phases);
            _problem          = problem;
            _constraintHandle = ConstraintFunction;
            Status            = OptimStatus.CREATED;
        }

        public int    N                { get; set; } = 6;
        public int    Maxits           { get; set; } = 4000;
        public double Epsx             { get; set; } = 0;
        public double Diffstep         { get; set; } = 1e-9;
        public double Stpmax           { get; set; } = 1;
        public int    OptimizerTimeout { get; set; } = 500000; // milliseconds

        public void Dispose()
        {
            // FIXME: ObjectPooling
        }

        private void CalculatePrimalFeasibility(double[] f)
        {
            Cost = f[0];
            // FIXME: skip inactive inequality constraints
            PrimalFeasibility = 0;
            for (int i = 1; i < f.Length; i++)
                PrimalFeasibility += f[i] * f[i];

            PrimalFeasibility = Sqrt(PrimalFeasibility);
        }

        private void ConstraintFunction(double[] x, double[] f, alglib.sparsematrix j, object? o)
        {
            _timeoutToken.ThrowIfCancellationRequested();

            var            lastPhase           = new VariableLayout(K, _numPhases - 1);
            VariableLayout optimizedPhase      = default;
            int            optimizedPhaseIndex = -1;

            for (int p = 0; p < _numPhases; p++)
            {
                if (_phases[p].Coast || !_phases[p].OptimizeTime)
                    continue;

                optimizedPhase      = new VariableLayout(K, p);
                optimizedPhaseIndex = p;
            }

            int ci = 0;

            // cost metric
            switch (_problem.Cost)
            {
                case Problem.ProblemCost.MIN_TIME:
                    if (optimizedPhaseIndex == -1)
                        throw new Exception("no phase with optimized time for MIN_TIME cost function");

                    f[ci++] = x[optimizedPhase.Tf] - x[optimizedPhase.Ti];
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, optimizedPhase.Ti, -1.0);
                    alglib.sparseappendelement(j, optimizedPhase.Tf, 1.0);

                    break;
                case Problem.ProblemCost.MAX_MASS:
                    if (optimizedPhaseIndex == -1)
                        throw new Exception("no phase with optimized time for MAX_MASS cost function");

                    f[ci++] = -x[optimizedPhase.MEnd];
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, optimizedPhase.MEnd, -1);

                    break;
                case Problem.ProblemCost.MAX_ENERGY:
                    var    vf   = V3.CopyFromIndices(x, lastPhase.VEnd);
                    var    rf   = V3.CopyFromIndices(x, lastPhase.REnd);
                    double rfm  = rf.magnitude;
                    double rfm3 = rfm * rfm * rfm;
                    // XXX: implicit mu = 1.0 on the next line
                    f[ci++] = -(0.5 * vf.sqrMagnitude - 1.0 / rfm);
                    alglib.sparseappendemptyrow(j);
                    // partials of 1/|r|
                    alglib.sparseappendelement(j, lastPhase.REnd.Item1, -1.0 / rfm3 * rf.x);
                    alglib.sparseappendelement(j, lastPhase.REnd.Item2, -1.0 / rfm3 * rf.y);
                    alglib.sparseappendelement(j, lastPhase.REnd.Item3, -1.0 / rfm3 * rf.z);
                    // partials of -1/2*|v|^2
                    alglib.sparseappendelement(j, lastPhase.VEnd.Item1, -vf.x);
                    alglib.sparseappendelement(j, lastPhase.VEnd.Item2, -vf.y);
                    alglib.sparseappendelement(j, lastPhase.VEnd.Item3, -vf.z);
                    break;
                default:
                    throw new Exception("code should not be reachable");
            }

            // initial conditions
            var firstPhase = new VariableLayout(K, 0);
            var r0         = V3.CopyFromIndices(x, firstPhase.RStart);
            var v0         = V3.CopyFromIndices(x, firstPhase.VStart);

            f[ci++] = r0.x - _problem.R0.x;
            f[ci++] = r0.y - _problem.R0.y;
            f[ci++] = r0.z - _problem.R0.z;
            f[ci++] = v0.x - _problem.V0.x;
            f[ci++] = v0.y - _problem.V0.y;
            f[ci++] = v0.z - _problem.V0.z;
            f[ci++] = x[firstPhase.Ti];

            alglib.sparseappendemptyrow(j);
            alglib.sparseappendelement(j, firstPhase.RStart.Item1, 1.0);
            alglib.sparseappendemptyrow(j);
            alglib.sparseappendelement(j, firstPhase.RStart.Item2, 1.0);
            alglib.sparseappendemptyrow(j);
            alglib.sparseappendelement(j, firstPhase.RStart.Item3, 1.0);
            alglib.sparseappendemptyrow(j);
            alglib.sparseappendelement(j, firstPhase.VStart.Item1, 1.0);
            alglib.sparseappendemptyrow(j);
            alglib.sparseappendelement(j, firstPhase.VStart.Item2, 1.0);
            alglib.sparseappendemptyrow(j);
            alglib.sparseappendelement(j, firstPhase.VStart.Item3, 1.0);
            alglib.sparseappendemptyrow(j);
            alglib.sparseappendelement(j, firstPhase.Ti, 1.0);

            for (int p = 0; p < _numPhases; p++)
            {
                var thisPhase = new VariableLayout(K, p);

                double mdot   = _phases[p].mdot;
                double thrust = _phases[p].thrust;

                double h = (x[thisPhase.Tf] - x[thisPhase.Ti]) / (N - 1);

                // dynamical constraints per phase
                for (int n = 0; n < N - 1; n += 1)
                {
                    (int rkx, int rky, int rkz) = thisPhase.R(2 * n);
                    (int vkx, int vky, int vkz) = thisPhase.V(2 * n);
                    int mk = thisPhase.M(2 * n);
                    (int ukx, int uky, int ukz) = thisPhase.U(2 * n);

                    // dr_x/dt = v_x
                    f[ci++] = x[rkx + 2] - x[rkx] - 1.0 / 6.0 * h * (x[vkx] + 4 * x[vkx + 1] + x[vkx + 2]);
                    f[ci++] = x[rkx + 1] - 0.5 * (x[rkx] + x[rkx + 2]) - 0.125 * h * (x[vkx] - x[vkx + 2]);
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rkx, -1.0);
                    alglib.sparseappendelement(j, rkx + 2, 1.0);
                    alglib.sparseappendelement(j, vkx, -1.0 / 6.0 * h);
                    alglib.sparseappendelement(j, vkx + 1, -4.0 / 6.0 * h);
                    alglib.sparseappendelement(j, vkx + 2, -1.0 / 6.0 * h);
                    alglib.sparseappendelement(j, thisPhase.Ti, 1.0 / 6.0 / (N - 1) * (x[vkx] + 4 * x[vkx + 1] + x[vkx + 2]));
                    alglib.sparseappendelement(j, thisPhase.Tf, -1.0 / 6.0 / (N - 1) * (x[vkx] + 4 * x[vkx + 1] + x[vkx + 2]));
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rkx, -0.5);
                    alglib.sparseappendelement(j, rkx + 1, 1.0);
                    alglib.sparseappendelement(j, rkx + 2, -0.5);
                    alglib.sparseappendelement(j, vkx, -0.125 * h);
                    alglib.sparseappendelement(j, vkx + 2, 0.125 * h);
                    alglib.sparseappendelement(j, thisPhase.Ti, 0.125 / (N - 1) * (x[vkx] - x[vkx + 2]));
                    alglib.sparseappendelement(j, thisPhase.Tf, -0.125 / (N - 1) * (x[vkx] - x[vkx + 2]));

                    // dr_y/dt = v_y
                    f[ci++] = x[rky + 2] - x[rky] - 1.0 / 6.0 * h * (x[vky] + 4 * x[vky + 1] + x[vky + 2]);
                    f[ci++] = x[rky + 1] - 0.5 * (x[rky] + x[rky + 2]) - 0.125 * h * (x[vky] - x[vky + 2]);
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rky, -1.0);
                    alglib.sparseappendelement(j, rky + 2, 1.0);
                    alglib.sparseappendelement(j, vky, -1.0 / 6.0 * h);
                    alglib.sparseappendelement(j, vky + 1, -4.0 / 6.0 * h);
                    alglib.sparseappendelement(j, vky + 2, -1.0 / 6.0 * h);
                    alglib.sparseappendelement(j, thisPhase.Ti, 1.0 / 6.0 / (N - 1) * (x[vky] + 4 * x[vky + 1] + x[vky + 2]));
                    alglib.sparseappendelement(j, thisPhase.Tf, -1.0 / 6.0 / (N - 1) * (x[vky] + 4 * x[vky + 1] + x[vky + 2]));
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rky, -0.5);
                    alglib.sparseappendelement(j, rky + 1, 1.0);
                    alglib.sparseappendelement(j, rky + 2, -0.5);
                    alglib.sparseappendelement(j, vky, -0.125 * h);
                    alglib.sparseappendelement(j, vky + 2, 0.125 * h);
                    alglib.sparseappendelement(j, thisPhase.Ti, 0.125 / (N - 1) * (x[vky] - x[vky + 2]));
                    alglib.sparseappendelement(j, thisPhase.Tf, -0.125 / (N - 1) * (x[vky] - x[vky + 2]));

                    // dr_z/dt = v_z
                    f[ci++] = x[rkz + 2] - x[rkz] - 1.0 / 6.0 * h * (x[vkz] + 4 * x[vkz + 1] + x[vkz + 2]);
                    f[ci++] = x[rkz + 1] - 0.5 * (x[rkz] + x[rkz + 2]) - 0.125 * h * (x[vkz] - x[vkz + 2]);
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rkz, -1.0);
                    alglib.sparseappendelement(j, rkz + 2, 1.0);
                    alglib.sparseappendelement(j, vkz, -1.0 / 6.0 * h);
                    alglib.sparseappendelement(j, vkz + 1, -4.0 / 6.0 * h);
                    alglib.sparseappendelement(j, vkz + 2, -1.0 / 6.0 * h);
                    alglib.sparseappendelement(j, thisPhase.Ti, 1.0 / 6.0 / (N - 1) * (x[vkz] + 4 * x[vkz + 1] + x[vkz + 2]));
                    alglib.sparseappendelement(j, thisPhase.Tf, -1.0 / 6.0 / (N - 1) * (x[vkz] + 4 * x[vkz + 1] + x[vkz + 2]));
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rkz, -0.5);
                    alglib.sparseappendelement(j, rkz + 1, 1.0);
                    alglib.sparseappendelement(j, rkz + 2, -0.5);
                    alglib.sparseappendelement(j, vkz, -0.125 * h);
                    alglib.sparseappendelement(j, vkz + 2, 0.125 * h);
                    alglib.sparseappendelement(j, thisPhase.Ti, 0.125 / (N - 1) * (x[vkz] - x[vkz + 2]));
                    alglib.sparseappendelement(j, thisPhase.Tf, -0.125 / (N - 1) * (x[vkz] - x[vkz + 2]));

                    double rm  = new V3(x[rkx], x[rky], x[rkz]).magnitude;
                    double r2  = rm * rm;
                    double r3  = r2 * rm;
                    double r5  = r3 * r2;
                    double r1m = new V3(x[rkx + 1], x[rky + 1], x[rkz + 1]).magnitude;
                    double r12 = r1m * r1m;
                    double r13 = r12 * r1m;
                    double r15 = r13 * r12;
                    double r2m = new V3(x[rkx + 2], x[rky + 2], x[rkz + 2]).magnitude;
                    double r22 = r2m * r2m;
                    double r23 = r22 * r2m;
                    double r25 = r23 * r22;

                    // dv_x/dt = -r_x/r3 + T/m*u_x
                    f[ci++] = x[vkx + 2] - x[vkx] - 1.0 / 6.0 * h * (-x[rkx] / r3 + thrust / x[mk] * x[ukx] + 4 * (-x[rkx + 1] / r13 + thrust / x[mk + 1] * x[ukx + 1]) - x[rkx + 2] / r23 + thrust / x[mk + 2] * x[ukx + 2]);
                    f[ci++] = x[vkx + 1] - 0.5 * (x[vkx] + x[vkx + 2]) - 0.125 * h * (-x[rkx] / r3 + thrust / x[mk] * x[ukx] - (-x[rkx + 2] / r23 + thrust / x[mk + 2] * x[ukx + 2]));
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rkx, -1.0 / 6.0 * h * (-1.0 / r3 + 3.0 * x[rkx] * x[rkx] / r5));
                    alglib.sparseappendelement(j, rkx + 1, -4.0 / 6.0 * h * (-1.0 / r13 + 3.0 * x[rkx + 1] * x[rkx + 1] / r15));
                    alglib.sparseappendelement(j, rkx + 2, -1.0 / 6.0 * h * (-1.0 / r23 + 3.0 * x[rkx + 2] * x[rkx + 2] / r25));
                    alglib.sparseappendelement(j, rky, -1.0 / 6.0 * h * (3.0 * x[rkx] * x[rky] / r5));
                    alglib.sparseappendelement(j, rky + 1, -4.0 / 6.0 * h * (3.0 * x[rkx + 1] * x[rky + 1] / r15));
                    alglib.sparseappendelement(j, rky + 2, -1.0 / 6.0 * h * (3.0 * x[rkx + 2] * x[rky + 2] / r25));
                    alglib.sparseappendelement(j, rkz, -1.0 / 6.0 * h * (3.0 * x[rkx] * x[rkz] / r5));
                    alglib.sparseappendelement(j, rkz + 1, -4.0 / 6.0 * h * (3.0 * x[rkx + 1] * x[rkz + 1] / r15));
                    alglib.sparseappendelement(j, rkz + 2, -1.0 / 6.0 * h * (3.0 * x[rkx + 2] * x[rkz + 2] / r25));
                    alglib.sparseappendelement(j, vkx, -1.0);
                    alglib.sparseappendelement(j, vkx + 2, 1.0);
                    alglib.sparseappendelement(j, mk, 1.0 / 6.0 * h * thrust * x[ukx] / (x[mk] * x[mk]));
                    alglib.sparseappendelement(j, mk + 1, 4.0 / 6.0 * h * thrust * x[ukx] / (x[mk + 1] * x[mk + 1]));
                    alglib.sparseappendelement(j, mk + 2, 1.0 / 6.0 * h * thrust * x[ukx] / (x[mk + 2] * x[mk + 2]));
                    alglib.sparseappendelement(j, ukx, -1.0 / 6.0 * h * thrust / x[mk]);
                    alglib.sparseappendelement(j, ukx + 1, -4.0 / 6.0 * h * thrust / x[mk + 1]);
                    alglib.sparseappendelement(j, ukx + 2, -1.0 / 6.0 * h * thrust / x[mk + 2]);
                    alglib.sparseappendelement(j, thisPhase.Ti, 1.0 / 6.0 / (N - 1) * (-x[rkx] / r3 + thrust / x[mk] * x[ukx] + 4 * (-x[rkx + 1] / r13 + thrust / x[mk + 1] * x[ukx + 1]) - x[rkx + 2] / r23 + thrust / x[mk + 2] * x[ukx + 2]));
                    alglib.sparseappendelement(j, thisPhase.Tf, -1.0 / 6.0 / (N - 1) * (-x[rkx] / r3 + thrust / x[mk] * x[ukx] + 4 * (-x[rkx + 1] / r13 + thrust / x[mk + 1] * x[ukx + 1]) - x[rkx + 2] / r23 + thrust / x[mk + 2] * x[ukx + 2]));
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rkx, -0.125 * h * (-1.0 / r3 + 3.0 * x[rkx] * x[rkx] / r5));
                    alglib.sparseappendelement(j, rkx + 2, 0.125 * h * (-1.0 / r23 + 3.0 * x[rkx + 2] * x[rkx + 2] / r25));
                    alglib.sparseappendelement(j, rky, -0.125 * h * (3.0 * x[rkx] * x[rky] / r5));
                    alglib.sparseappendelement(j, rky + 2, 0.125 * h * (3.0 * x[rkx + 2] * x[rky + 2] / r25));
                    alglib.sparseappendelement(j, rkz, -0.125 * h * (3.0 * x[rkx] * x[rkz] / r5));
                    alglib.sparseappendelement(j, rkz + 2, 0.125 * h * (3.0 * x[rkx + 2] * x[rkz + 2] / r25));
                    alglib.sparseappendelement(j, vkx, -0.5);
                    alglib.sparseappendelement(j, vkx + 1, 1.0);
                    alglib.sparseappendelement(j, vkx + 2, -0.5);
                    alglib.sparseappendelement(j, mk, 0.125 * h * thrust * x[ukx] / (x[mk] * x[mk]));
                    alglib.sparseappendelement(j, mk + 2, -0.125 * h * thrust * x[ukx + 2] / (x[mk + 2] * x[mk + 2]));
                    alglib.sparseappendelement(j, ukx, -0.125 * h * thrust / x[mk]);
                    alglib.sparseappendelement(j, ukx + 2, 0.125 * h * thrust / x[mk + 2]);
                    alglib.sparseappendelement(j, thisPhase.Ti, 0.125 / (N - 1) * (-x[rkx] / r3 + thrust / x[mk] * x[ukx] - (-x[rkx + 2] / r23 + thrust / x[mk + 2] * x[ukx + 2])));
                    alglib.sparseappendelement(j, thisPhase.Tf, -0.125 / (N - 1) * (-x[rkx] / r3 + thrust / x[mk] * x[ukx] - (-x[rkx + 2] / r23 + thrust / x[mk + 2] * x[ukx + 2])));

                    // dv_y/dt = -r_y/r3 + T/m*u_y
                    f[ci++] = x[vky + 2] - x[vky] - 1.0 / 6.0 * h * (-x[rky] / r3 + thrust / x[mk] * x[uky] + 4 * (-x[rky + 1] / r13 + thrust / x[mk + 1] * x[uky + 1]) - x[rky + 2] / r23 + thrust / x[mk + 2] * x[uky + 2]);
                    f[ci++] = x[vky + 1] - 0.5 * (x[vky] + x[vky + 2]) - 0.125 * h * (-x[rky] / r3 + thrust / x[mk] * x[uky] - (-x[rky + 2] / r23 + thrust / x[mk + 2] * x[uky + 2]));
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rkx, -1.0 / 6.0 * h * (3.0 * x[rky] * x[rkx] / r5));
                    alglib.sparseappendelement(j, rkx + 1, -4.0 / 6.0 * h * (3.0 * x[rky + 1] * x[rkx + 1] / r15));
                    alglib.sparseappendelement(j, rkx + 2, -1.0 / 6.0 * h * (3.0 * x[rky + 2] * x[rkx + 2] / r25));
                    alglib.sparseappendelement(j, rky, -1.0 / 6.0 * h * (-1.0 / r3 + 3.0 * x[rky] * x[rky] / r5));
                    alglib.sparseappendelement(j, rky + 1, -4.0 / 6.0 * h * (-1.0 / r13 + 3.0 * x[rky + 1] * x[rky + 1] / r15));
                    alglib.sparseappendelement(j, rky + 2, -1.0 / 6.0 * h * (-1.0 / r23 + 3.0 * x[rky + 2] * x[rky + 2] / r25));
                    alglib.sparseappendelement(j, rkz, -1.0 / 6.0 * h * (3.0 * x[rky] * x[rkz] / r5));
                    alglib.sparseappendelement(j, rkz + 1, -4.0 / 6.0 * h * (3.0 * x[rky + 1] * x[rkz + 1] / r15));
                    alglib.sparseappendelement(j, rkz + 2, -1.0 / 6.0 * h * (3.0 * x[rky + 2] * x[rkz + 2] / r25));
                    alglib.sparseappendelement(j, vky, -1.0);
                    alglib.sparseappendelement(j, vky + 2, 1.0);
                    alglib.sparseappendelement(j, mk, 1.0 / 6.0 * h * thrust * x[uky] / (x[mk] * x[mk]));
                    alglib.sparseappendelement(j, mk + 1, 4.0 / 6.0 * h * thrust * x[uky + 1] / (x[mk + 1] * x[mk + 1]));
                    alglib.sparseappendelement(j, mk + 2, 1.0 / 6.0 * h * thrust * x[uky + 2] / (x[mk + 2] * x[mk + 2]));
                    alglib.sparseappendelement(j, uky, -1.0 / 6.0 * h * thrust / x[mk]);
                    alglib.sparseappendelement(j, uky + 1, -4.0 / 6.0 * h * thrust / x[mk + 1]);
                    alglib.sparseappendelement(j, uky + 2, -1.0 / 6.0 * h * thrust / x[mk + 2]);
                    alglib.sparseappendelement(j, thisPhase.Ti, 1.0 / 6.0 / (N - 1) * (-x[rky] / r3 + thrust / x[mk] * x[uky] + 4 * (-x[rky + 1] / r13 + thrust / x[mk + 1] * x[uky + 1]) - x[rky + 2] / r23 + thrust / x[mk + 2] * x[uky + 2]));
                    alglib.sparseappendelement(j, thisPhase.Tf, -1.0 / 6.0 / (N - 1) * (-x[rky] / r3 + thrust / x[mk] * x[uky] + 4 * (-x[rky + 1] / r13 + thrust / x[mk + 1] * x[uky + 1]) - x[rky + 2] / r23 + thrust / x[mk + 2] * x[uky + 2]));
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rkx, -0.125 * h * (3.0 * x[rky] * x[rkx] / r5));
                    alglib.sparseappendelement(j, rkx + 2, 0.125 * h * (3.0 * x[rky + 2] * x[rkx + 2] / r25));
                    alglib.sparseappendelement(j, rky, -0.125 * h * (-1.0 / r3 + 3.0 * x[rky] * x[rky] / r5));
                    alglib.sparseappendelement(j, rky + 2, 0.125 * h * (-1.0 / r23 + 3.0 * x[rky + 2] * x[rky + 2] / r25));
                    alglib.sparseappendelement(j, rkz, -0.125 * h * (3.0 * x[rky] * x[rkz] / r5));
                    alglib.sparseappendelement(j, rkz + 2, 0.125 * h * (3.0 * x[rky + 2] * x[rkz + 2] / r25));
                    alglib.sparseappendelement(j, vky, -0.5);
                    alglib.sparseappendelement(j, vky + 1, 1.0);
                    alglib.sparseappendelement(j, vky + 2, -0.5);
                    alglib.sparseappendelement(j, mk, 0.125 * h * thrust * x[uky] / (x[mk] * x[mk]));
                    alglib.sparseappendelement(j, mk + 2, -0.125 * h * thrust * x[uky + 2] / (x[mk + 2] * x[mk + 2]));
                    alglib.sparseappendelement(j, uky, -0.125 * h * thrust / x[mk]);
                    alglib.sparseappendelement(j, uky + 2, 0.125 * h * thrust / x[mk + 2]);
                    alglib.sparseappendelement(j, thisPhase.Ti, 0.125 / (N - 1) * (-x[rky] / r3 + thrust / x[mk] * x[uky] - (-x[rky + 2] / r23 + thrust / x[mk + 2] * x[uky + 2])));
                    alglib.sparseappendelement(j, thisPhase.Tf, -0.125 / (N - 1) * (-x[rky] / r3 + thrust / x[mk] * x[uky] - (-x[rky + 2] / r23 + thrust / x[mk + 2] * x[uky + 2])));

                    // dv_z/dt = -r_z/r3 + T/m*u_z
                    f[ci++] = x[vkz + 2] - x[vkz] - 1.0 / 6.0 * h * (-x[rkz] / r3 + thrust / x[mk] * x[ukz] + 4 * (-x[rkz + 1] / r13 + thrust / x[mk + 1] * x[ukz + 1]) - x[rkz + 2] / r23 + thrust / x[mk + 2] * x[ukz + 2]);
                    f[ci++] = x[vkz + 1] - 0.5 * (x[vkz] + x[vkz + 2]) - 0.125 * h * (-x[rkz] / r3 + thrust / x[mk] * x[ukz] - (-x[rkz + 2] / r23 + thrust / x[mk + 2] * x[ukz + 2]));
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rkx, -1.0 / 6.0 * h * (3.0 * x[rkz] * x[rkx] / r5));
                    alglib.sparseappendelement(j, rkx + 1, -4.0 / 6.0 * h * (3.0 * x[rkz + 1] * x[rkx + 1] / r15));
                    alglib.sparseappendelement(j, rkx + 2, -1.0 / 6.0 * h * (3.0 * x[rkz + 2] * x[rkx + 2] / r25));
                    alglib.sparseappendelement(j, rky, -1.0 / 6.0 * h * (3.0 * x[rkz] * x[rky] / r5));
                    alglib.sparseappendelement(j, rky + 1, -4.0 / 6.0 * h * (3.0 * x[rkz + 1] * x[rky + 1] / r15));
                    alglib.sparseappendelement(j, rky + 2, -1.0 / 6.0 * h * (3.0 * x[rkz + 2] * x[rky + 2] / r25));
                    alglib.sparseappendelement(j, rkz, -1.0 / 6.0 * h * (-1.0 / r3 + 3.0 * x[rkz] * x[rkz] / r5));
                    alglib.sparseappendelement(j, rkz + 1, -4.0 / 6.0 * h * (-1.0 / r13 + 3.0 * x[rkz + 1] * x[rkz + 1] / r15));
                    alglib.sparseappendelement(j, rkz + 2, -1.0 / 6.0 * h * (-1.0 / r23 + 3.0 * x[rkz + 2] * x[rkz + 2] / r25));
                    alglib.sparseappendelement(j, vkz, -1.0);
                    alglib.sparseappendelement(j, vkz + 2, 1.0);
                    alglib.sparseappendelement(j, mk, 1.0 / 6.0 * h * thrust * x[ukz] / (x[mk] * x[mk]));
                    alglib.sparseappendelement(j, mk + 1, 4.0 / 6.0 * h * thrust * x[ukz + 1] / (x[mk + 1] * x[mk + 1]));
                    alglib.sparseappendelement(j, mk + 2, 1.0 / 6.0 * h * thrust * x[ukz + 2] / (x[mk + 2] * x[mk + 2]));
                    alglib.sparseappendelement(j, ukz, -1.0 / 6.0 * h * thrust / x[mk]);
                    alglib.sparseappendelement(j, ukz + 1, -4.0 / 6.0 * h * thrust / x[mk + 1]);
                    alglib.sparseappendelement(j, ukz + 2, -1.0 / 6.0 * h * thrust / x[mk + 2]);
                    alglib.sparseappendelement(j, thisPhase.Ti, 1.0 / 6.0 / (N - 1) * (-x[rkz] / r3 + thrust / x[mk] * x[ukz] + 4 * (-x[rkz + 1] / r13 + thrust / x[mk + 1] * x[ukz + 1]) - x[rkz + 2] / r23 + thrust / x[mk + 2] * x[ukz + 2]));
                    alglib.sparseappendelement(j, thisPhase.Tf, -1.0 / 6.0 / (N - 1) * (-x[rkz] / r3 + thrust / x[mk] * x[ukz] + 4 * (-x[rkz + 1] / r13 + thrust / x[mk + 1] * x[ukz + 1]) - x[rkz + 2] / r23 + thrust / x[mk + 2] * x[ukz + 2]));
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, rkx, -0.125 * h * (3.0 * x[rkz] * x[rkx] / r5));
                    alglib.sparseappendelement(j, rkx + 2, 0.125 * h * (3.0 * x[rkz + 2] * x[rkx + 2] / r25));
                    alglib.sparseappendelement(j, rky, -0.125 * h * (3.0 * x[rkz] * x[rky] / r5));
                    alglib.sparseappendelement(j, rky + 2, 0.125 * h * (3.0 * x[rkz + 2] * x[rky + 2] / r25));
                    alglib.sparseappendelement(j, rkz, -0.125 * h * (-1.0 / r3 + 3.0 * x[rkz] * x[rkz] / r5));
                    alglib.sparseappendelement(j, rkz + 2, 0.125 * h * (-1.0 / r23 + 3.0 * x[rkz + 2] * x[rkz + 2] / r25));
                    alglib.sparseappendelement(j, vkz, -0.5);
                    alglib.sparseappendelement(j, vkz + 1, 1.0);
                    alglib.sparseappendelement(j, vkz + 2, -0.5);
                    alglib.sparseappendelement(j, mk, 0.125 * h * thrust * x[ukz] / (x[mk] * x[mk]));
                    alglib.sparseappendelement(j, mk + 2, -0.125 * h * thrust * x[ukz + 2] / (x[mk + 2] * x[mk + 2]));
                    alglib.sparseappendelement(j, ukz, -0.125 * h * thrust / x[mk]);
                    alglib.sparseappendelement(j, ukz + 2, 0.125 * h * thrust / x[mk + 2]);
                    alglib.sparseappendelement(j, thisPhase.Ti, 0.125 / (N - 1) * (-x[rkz] / r3 + thrust / x[mk] * x[ukz] - (-x[rkz + 2] / r23 + thrust / x[mk + 2] * x[ukz + 2])));
                    alglib.sparseappendelement(j, thisPhase.Tf, -0.125 / (N - 1) * (-x[rkz] / r3 + thrust / x[mk] * x[ukz] - (-x[rkz + 2] / r23 + thrust / x[mk + 2] * x[ukz + 2])));

                    // dm/dt = -mdot
                    f[ci++] = x[mk + 2] - x[mk] + h * mdot;
                    f[ci++] = x[mk + 1] - 0.5 * (x[mk] + x[mk + 2]);
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, mk, -1.0);
                    alglib.sparseappendelement(j, mk + 2, 1.0);
                    alglib.sparseappendelement(j, thisPhase.Ti, -1.0 / (N - 1) * mdot);
                    alglib.sparseappendelement(j, thisPhase.Tf, 1.0 / (N - 1) * mdot);
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, mk, -0.5);
                    alglib.sparseappendelement(j, mk + 1, 1.0);
                    alglib.sparseappendelement(j, mk + 2, -0.5);
                }

                // mass constraints per phase
                f[ci++] = x[thisPhase.MStart] - _phases[p].m0;
                alglib.sparseappendemptyrow(j);
                alglib.sparseappendelement(j, thisPhase.MStart, 1.0);

                // burn time constraints per phase
                if (_phases[p].OptimizeTime)
                {
                    f[ci++] = 0;
                    alglib.sparseappendemptyrow(j);
                }
                else
                {
                    f[ci++] = x[thisPhase.Tf] - x[thisPhase.Ti] - _phases[p].bt;
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, thisPhase.Ti, -1.0);
                    alglib.sparseappendelement(j, thisPhase.Tf, 1.0);
                }

                // path inequality control constraints
                for (int k = 0; k < K; k++)
                {
                    var u = V3.CopyFromIndices(x, thisPhase.U(k));
                    f[ci++] = u.sqrMagnitude - 1.0;
                    alglib.sparseappendemptyrow(j);
                    alglib.sparseappendelement(j, thisPhase.U(k).Item1, 2 * u.x);
                    alglib.sparseappendelement(j, thisPhase.U(k).Item2, 2 * u.y);
                    alglib.sparseappendelement(j, thisPhase.U(k).Item3, 2 * u.z);
                }

                // path equality control constraints (for unguided phases)
                if (_phases[p].Unguided)
                {
                    if (p == 0)
                    {
                        f[ci++] = 0;
                        f[ci++] = 0;
                        f[ci++] = 0;
                        alglib.sparseappendemptyrow(j);
                        alglib.sparseappendemptyrow(j);
                        alglib.sparseappendemptyrow(j);
                    }
                    else
                    {
                        var prevPhase = new VariableLayout(K, p - 1);

                        var u  = V3.CopyFromIndices(x, thisPhase.UStart);
                        var u0 = V3.CopyFromIndices(x, prevPhase.UEnd);

                        f[ci++] = u.x - u0.x;
                        f[ci++] = u.y - u0.y;
                        f[ci++] = u.z - u0.z;
                        alglib.sparseappendemptyrow(j);
                        alglib.sparseappendelement(j, prevPhase.UEnd.Item1, -1.0);
                        alglib.sparseappendelement(j, thisPhase.UStart.Item1, 1.0);
                        alglib.sparseappendemptyrow(j);
                        alglib.sparseappendelement(j, prevPhase.UEnd.Item2, -1.0);
                        alglib.sparseappendelement(j, thisPhase.UStart.Item2, 1.0);
                        alglib.sparseappendemptyrow(j);
                        alglib.sparseappendelement(j, prevPhase.UEnd.Item3, -1.0);
                        alglib.sparseappendelement(j, thisPhase.UStart.Item3, 1.0);
                    }

                    for (int k = 1; k < K; k++)
                    {
                        var u     = V3.CopyFromIndices(x, thisPhase.U(k));
                        var uPrev = V3.CopyFromIndices(x, thisPhase.U(k - 1));
                        f[ci++] = u.x - uPrev.x;
                        f[ci++] = u.y - uPrev.y;
                        f[ci++] = u.z - uPrev.z;
                        alglib.sparseappendemptyrow(j);
                        alglib.sparseappendelement(j, thisPhase.U(k - 1).Item1, -1.0);
                        alglib.sparseappendelement(j, thisPhase.U(k).Item1, 1.0);
                        alglib.sparseappendemptyrow(j);
                        alglib.sparseappendelement(j, thisPhase.U(k - 1).Item2, -1.0);
                        alglib.sparseappendelement(j, thisPhase.U(k).Item2, 1.0);
                        alglib.sparseappendemptyrow(j);
                        alglib.sparseappendelement(j, thisPhase.U(k - 1).Item3, -1.0);
                        alglib.sparseappendelement(j, thisPhase.U(k).Item3, 1.0);
                    }
                }
                else
                {
                    for (int k = 0; k < K; k++)
                    {
                        f[ci++] = 0;
                        f[ci++] = 0;
                        f[ci++] = 0;
                        alglib.sparseappendemptyrow(j);
                        alglib.sparseappendemptyrow(j);
                        alglib.sparseappendemptyrow(j);
                    }
                }

                if (p == _numPhases - 1) break;

                // continuity conditions per phase
                var nextPhase = new VariableLayout(K, p + 1);

                var rMinus = V3.CopyFromIndices(x, thisPhase.REnd);
                var vMinus = V3.CopyFromIndices(x, thisPhase.VEnd);
                var rPlus  = V3.CopyFromIndices(x, nextPhase.RStart);
                var vPlus  = V3.CopyFromIndices(x, nextPhase.VStart);

                f[ci++] = rPlus.x - rMinus.x;
                f[ci++] = rPlus.y - rMinus.y;
                f[ci++] = rPlus.z - rMinus.z;
                f[ci++] = vPlus.x - vMinus.x;
                f[ci++] = vPlus.y - vMinus.y;
                f[ci++] = vPlus.z - vMinus.z;
                f[ci++] = x[nextPhase.Ti] - x[thisPhase.Tf];

                alglib.sparseappendemptyrow(j);
                alglib.sparseappendelement(j, thisPhase.REnd.Item1, -1.0);
                alglib.sparseappendelement(j, nextPhase.RStart.Item1, 1.0);
                alglib.sparseappendemptyrow(j);
                alglib.sparseappendelement(j, thisPhase.REnd.Item2, -1.0);
                alglib.sparseappendelement(j, nextPhase.RStart.Item2, 1.0);
                alglib.sparseappendemptyrow(j);
                alglib.sparseappendelement(j, thisPhase.REnd.Item3, -1.0);
                alglib.sparseappendelement(j, nextPhase.RStart.Item3, 1.0);
                alglib.sparseappendemptyrow(j);
                alglib.sparseappendelement(j, thisPhase.VEnd.Item1, -1.0);
                alglib.sparseappendelement(j, nextPhase.VStart.Item1, 1.0);
                alglib.sparseappendemptyrow(j);
                alglib.sparseappendelement(j, thisPhase.VEnd.Item2, -1.0);
                alglib.sparseappendelement(j, nextPhase.VStart.Item2, 1.0);
                alglib.sparseappendemptyrow(j);
                alglib.sparseappendelement(j, thisPhase.VEnd.Item3, -1.0);
                alglib.sparseappendelement(j, nextPhase.VStart.Item3, 1.0);
                alglib.sparseappendemptyrow(j);
                alglib.sparseappendelement(j, thisPhase.Tf, -1.0);
                alglib.sparseappendelement(j, nextPhase.Ti, 1.0);
            }

            // terminal constraints
            _problem.Terminal.Constraints(x, lastPhase.REnd, lastPhase.VEnd, f, j, ref ci);

            if (ci != VariableLayout.NumConstraintsForPhases(N, K, _numPhases, _problem) + 1)
                throw new Exception("Constraint num mismatch");
        }

        private class VacuumThrustKernel
        {
            public static int N => InterpolantLayout.INTERPOLANT_LAYOUT_LEN;

            public Phase Phase = null!;

            public void dydt(IList<double> yin, double x, IList<double> dyout)
            {
                Check.True(Phase.Normalized);

                var y  = InterpolantLayout.CreateFrom(yin);
                var dy = new InterpolantLayout();

                double at = Phase.thrust / y.M;

                double r2 = y.R.sqrMagnitude;
                double r  = Sqrt(r2);
                double r3 = r2 * r;

                dy.R  = y.V;
                dy.V  = -y.R / r3 + at * y.U;
                dy.M  = -Phase.mdot;
                dy.U  = V3.zero;
                dy.DV = at;

                dy.CopyTo(dyout);
            }
        }

        private readonly VacuumThrustKernel _ode    = new VacuumThrustKernel();
        private readonly DP5                _solver = new DP5();

        public void Integrate(Vn y0, Vn yf, Phase phase, double t0, double tf, Solution solution)
        {
            _solver.ThrowOnMaxIter = true;
            _solver.Maxiter        = 2000;
            _solver.Rtol           = 1e-6;
            _solver.Atol           = 1e-6;
            _ode.Phase             = phase;
            var interpolant = Hn.Get(VacuumThrustKernel.N);
            _solver.Solve(_ode.dydt, y0, yf, t0, tf, interpolant);
            solution.AddSegment(t0, tf, interpolant, phase);
        }

        private void Shooting(Solution solution, V3 u0)
        {
            using var initial  = Vn.Rent(InterpolantLayout.INTERPOLANT_LAYOUT_LEN);
            using var terminal = Vn.Rent(InterpolantLayout.INTERPOLANT_LAYOUT_LEN);

            var y0 = new InterpolantLayout();

            double t0 = 0;

            y0.R  = _problem.R0;
            y0.V  = _problem.V0;
            y0.U  = u0;
            y0.DV = 0;

            for (int p = 0; p < _numPhases; p++)
            {
                Phase phase = _phases[p];
                y0.M = phase.m0;

                double bt = phase.bt;
                double tf = t0 + bt;

                y0.CopyTo(initial);

                Integrate(initial, terminal, phase, t0, tf, solution);

                y0.CopyFrom(terminal);

                t0 += bt;
            }
        }

        // XXX: tightly coupled to an ascent
        // but needs the Problem and the normalized/scaled variables, which are more accessible here for now
        public Solution InitialGuess(double incT)
        {
            V3 r0 = _problem.R0;
            // guess the initial launch direction
            V3 enu = Astro.ENUHeadingForInclination(incT, r0);
            // add 45 degrees up
            enu.z = 1.0;
            V3 u0 = Astro.ENUToECI(r0, enu).normalized;

            var solution = new Solution(_problem);

            Shooting(solution, u0);

            return solution;
        }

        private void TranscribeGuess(double[] x, Solution oldSolution)
        {
            double tOffset = 0; // TODO: time offset into old solution
            double t0      = 0;
            for (int p = 0; p < _numPhases; p++)
            {
                // WARNING: phases in the old solution may not match phases in the current problem
                Phase phase     = _phases[p];
                var   thisPhase = new VariableLayout(K, p);

                // caller needs to ensure burn times in the phases are accurate -- see Ascent.ApplyOldBurnTimesToPhases()
                double tf = t0 + phase.bt;
                double h  = (tf - t0) / (K - 1);
                double m0 = phase.m0;

                for (int k = 0; k < K; k++)
                {
                    double bt = k * h;
                    double t  = t0 + bt;

                    V3 r = oldSolution.RBar(t + tOffset);
                    r.CopyToIndices(x, thisPhase.R(k));

                    V3 v = oldSolution.VBar(t + tOffset);
                    v.CopyToIndices(x, thisPhase.V(k));

                    V3 u = oldSolution.UBar(t + tOffset);
                    u.CopyToIndices(x, thisPhase.U(k));

                    x[thisPhase.M(k)] = m0 - phase.mdot * bt;
                }

                x[thisPhase.Ti] = t0;
                x[thisPhase.Tf] = tf;

                t0 += phase.bt;
            }
        }

        private Solution UnSafeRun(Solution oldSolution)
        {
            double[] xGuess = new double[VariableLayout.NumVariablesForPhases(K, _numPhases)];
            double[] x      = new double[VariableLayout.NumVariablesForPhases(K, _numPhases)];
            double[] bndl   = new double[VariableLayout.NumVariablesForPhases(K, _numPhases)];
            double[] bndu   = new double[VariableLayout.NumVariablesForPhases(K, _numPhases)];
            double[] f      = new double[VariableLayout.NumConstraintsForPhases(N, K, _numPhases, _problem) + 1];

            TranscribeGuess(xGuess, oldSolution);

            // XXX: debugging initial constraints, delete this
            //alglib.sparsecreatecrsempty(VariableLayout.NumVariablesForPhases(K, _numPhases), out alglib.sparsematrix j2);
            //ConstraintFunction(xGuess, f, j2, null);

            for (int i = 0; i < bndu.Length; i++)
            {
                bndu[i] = double.PositiveInfinity;
                bndl[i] = double.NegativeInfinity;
            }

            // box constraints on initial conditions

            var firstPhase = new VariableLayout(K, 0);
            (int r0X, int r0Y, int r0Z) = firstPhase.RStart;
            (int v0X, int v0Y, int v0Z) = firstPhase.VStart;

            bndu[r0X] = bndl[r0X] = _problem.R0.x;
            bndu[r0Y] = bndl[r0Y] = _problem.R0.y;
            bndu[r0Z] = bndl[r0Z] = _problem.R0.z;
            bndu[v0X] = bndl[v0X] = _problem.V0.x;
            bndu[v0Y] = bndl[v0Y] = _problem.V0.y;
            bndu[v0Z] = bndl[v0Z] = _problem.V0.z;

            // FIXME: set boundaries on fixed times up until we hit an optimized stage
            // FIXME: set box path boundaries on mass
            // FIXME: set box path boundaries on control

            for (int p = 0; p < _numPhases; p++)
            {
                // box constraints on initial stage masses
                var thisPhase = new VariableLayout(K, p);

                bndu[thisPhase.MStart] = bndl[thisPhase.MStart] = _phases[p].m0;
            }

            alglib.minnlccreate(VariableLayout.NumVariablesForPhases(K, _numPhases), xGuess, out _state);
            alglib.minnlcsetbc(_state, bndl, bndu);
            alglib.minnlcsetstpmax(_state, Stpmax);
            alglib.minnlcsetalgosqp(_state);
            alglib.minnlcsetcond(_state, Epsx, Maxits);
            alglib.minnlcsetnlc(_state, VariableLayout.NumConstraintsForPhases(N, K, _numPhases, _problem), 0);

            //alglib.minnlcoptguardgradient(_state, Diffstep);
            //alglib.minnlcoptguardsmoothness(_state, 1);

            alglib.minnlcoptimize(_state, _constraintHandle, null, null);
            alglib.minnlcresultsbuf(_state, ref x, _rep);

            TerminationType = _rep.terminationtype;
            Iterations      = _rep.iterationscount;

            DebugPrint("terminationtype: " + TerminationType);
            DebugPrint("iterations: " + Iterations);

            alglib.minnlcoptguardresults(_state, out alglib.optguardreport ogrep);

            if (ogrep.badgradsuspected)
                if (!DoubleMatrixSparsityValidation(ogrep.badgraduser, ogrep.badgradnum, 1e-3))
                    throw new Exception(
                        $"badgradsuspected: {ogrep.badgradfidx},{ogrep.badgradvidx}\nuser:\n{DoubleMatrixString(ogrep.badgraduser)}\nnumerical:\n{DoubleMatrixString(ogrep.badgradnum)}\nsparsity check:\n{DoubleMatrixSparsityCheck(ogrep.badgraduser, ogrep.badgradnum, 1e-3)}");

            if (ogrep.nonc0suspected)
                throw new Exception("nonc0suspected");

            if (ogrep.nonc1suspected)
                throw new Exception("nonc1suspected");

            alglib.sparsecreatecrsempty(VariableLayout.NumVariablesForPhases(K, _numPhases), out alglib.sparsematrix j);

            if (_rep.terminationtype != 8)
                ConstraintFunction(x, f, j, null);

            CalculatePrimalFeasibility(f);

            DebugPrint($"Cost: {Cost}");
            DebugPrint($"PrimalFeasibility: {PrimalFeasibility}");

            Status = PrimalFeasibility > 1e-6 ? OptimStatus.FAILED : OptimStatus.SUCCESS;

            return GenerateSolution(x);
        }

        public Solution? Run(Solution previousSolution)
        {
            foreach (Phase phase in _phases)
                DebugPrint(phase.ToString());

            try
            {
                var tokenSource = new CancellationTokenSource();
                tokenSource.CancelAfter(OptimizerTimeout);
                _timeoutToken = tokenSource.Token;
                Solution      = UnSafeRun(previousSolution);
                return Solution;
            }
            catch (OperationCanceledException)
            {
                Status = OptimStatus.CANCELLED;
            }

            return null;
        }

        private Solution GenerateSolution(double[] x)
        {
            var solution = new Solution(_problem);

            var temp  = new InterpolantLayout();                           // FIXME: rename
            var temp2 = Vn.Rent(InterpolantLayout.INTERPOLANT_LAYOUT_LEN); // FIXME: rename

            double dv = 0;

            for (int p = 0; p < _numPhases; p++)
            {
                var thisPhase   = new VariableLayout(K, p);
                var interpolant = Hn.Get(InterpolantLayout.INTERPOLANT_LAYOUT_LEN);

                double ti = x[thisPhase.Ti];
                double tf = x[thisPhase.Tf];
                double h  = (tf - ti) / (N - 1);

                for (int n = 0; n < N; n++)
                {
                    temp.R = V3.CopyFromIndices(x, thisPhase.R(2 * n));
                    temp.V = V3.CopyFromIndices(x, thisPhase.V(2 * n));
                    temp.M = x[thisPhase.M(2 * n)];
                    temp.U = V3.CopyFromIndices(x, thisPhase.U(2 * n));

                    double t = ti + n * h;
                    temp.DV = dv + _phases[p].DeltaVForTime(t - ti);

                    temp.CopyTo(temp2);
                    interpolant.Add(t, temp2);
                }

                solution.AddSegment(ti, tf, interpolant, _phases[p]);
                dv = temp.DV;
            }

            return solution;
        }

        public bool Success() => Status == OptimStatus.SUCCESS;

        public static OptimizerBuilder Builder() => new OptimizerBuilder();
    }
}
