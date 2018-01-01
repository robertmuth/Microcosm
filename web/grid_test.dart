import 'dart:math' as Math;
import 'dart:isolate';
import 'dart:collection';

import 'package:vector_math/vector_math.dart' as VM;

import 'grid.dart';
import 'shapes.dart' as Shapes;

const double epsilon = 0.000000001;

List<Shapes.Form> allShapes = [
  new Shapes.KnotAndSpheres(3, 2, 3),
  new Shapes.KnotAndSpheres(2, 3, 3),
  new Shapes.KnotAndSpheres(2, 5, 3),
  new Shapes.KnotAndSpheres(3, 4, 2),
  new Shapes.KnotAndSpheres(3, 5, 2),
  new Shapes.KnotAndSpheres(3, 4, 2),
  //
  new Shapes.HeartForm(),
  new Shapes.Metaballs(1),
  new Shapes.Metaballs(7),
  new Shapes.TriangleOfSpheres(8),
  new Shapes.StringOfEllipsoids(8, 0.7),
  new Shapes.Brain(3),
  new Shapes.Flower(14, 1.0),
  // KnotAndSpheres // not started
  // KnotAndTorus // not started
  new Shapes.Orbit(),
  new Shapes.TorusBox(),
  new Shapes.RingOfTori(5),
  new Shapes.SpheresAndCapsules(5),
  new Shapes.CubesAndCapsules(6, 0.6),
  // UFO
  new Shapes.Rings(5),
  new Shapes.Grinder(6),
  new Shapes.Kube(),
  new Shapes.Kube2(),
  new Shapes.Kube3(),
  new Shapes.Kube4(),
  new Shapes.Hexahedron(),
  new Shapes.Cylinder(),
  new Shapes.Octahedron(),
  new Shapes.Tetrahedron(),
  new Shapes.UFO(10),
  // Tennis
];

final bool reportSurfaceStats = true;
final bool compareOutput = false;

void DumpGridStats(
    Grid grid, DfsOption opt, int evalsBefore, MeshInfo info, int evalsAfter) {
  print("Evals ${evalsBefore}");
  List<int> outerSurface = grid.GetSurfacePoints(true, opt);
  print("Outer points: ${outerSurface.length}");
  //for (Point p in outerSurface) print("outer $p");

  List<int> innerSurface = grid.GetSurfacePoints(false, opt);
  print("Inner points: ${innerSurface.length}");
  //for (Point p in innerSurface) print("inner $p");
  print("Evals ${evalsAfter}");
  print("${info}");
}

MeshInfo AlgoBruteForce(
    Grid grid, DfsOption opt, Shapes.Form form, bool computeNormals) {
  form.evals = 0;
  grid.ValueUpdateBruteForce(opt);
  int evalsBefore = form.evals;
  List<int> hot = grid.ValidCubePoints();
  MeshInfo info = grid.CubeMarch(hot, opt, computeNormals);
  if (reportSurfaceStats) {
    print("BruteForce: (hot: ${hot.length})");
    DumpGridStats(grid, opt, evalsBefore, info, form.evals);
  }
  return info;
}

MeshInfo AlgoDfs(
    Grid grid, DfsOption opt, Shapes.Form form, bool computeNormals) {
  form.evals = 0;
  List<int> seeds = [];
  for (VM.Vector3 v in form.InteriorPoints()) {
    int s = grid.FindSurfacePointFromSeedBits(v, opt.eval);
    assert(s != -1);
    seeds.add(s);
  }
  List<int> hot = [];
  grid.ValueUpdateDfs(seeds, true, hot, opt);
  HashSet<int> hots = new HashSet.from(hot);
  assert(hots.length == hot.length);

  final int evalsBefore = form.evals;
  MeshInfo info = grid.CubeMarch(hot, opt, computeNormals);
  if (reportSurfaceStats) {
    print("Dfs: (hot: ${hot.length})");
    DumpGridStats(grid, opt, evalsBefore, info, form.evals);
  }
  return info;
}

MeshInfo AlgoRecursive(
    Grid grid, DfsOption opt, Shapes.Form form, bool computeNormals) {
  form.evals = 0;
  List<int> seeds = [];
  for (VM.Vector3 v in form.InteriorPoints()) {
    int s = grid.FindSurfacePointFromSeedBits(v, opt.eval);
    assert(s != -1, "point is ${v}");
    seeds.add(s);
  }
  List<int> hot = [];
  grid.ValueUpdateRecursive(seeds, opt.ttl, true, hot, opt);
  HashSet<int> hots = new HashSet.from(hot);
  assert(hots.length == hot.length);
  int evalsBefore = form.evals;
  MeshInfo info = grid.CubeMarch(hot, opt, computeNormals);
  if (reportSurfaceStats) {
    print("RecursiveDfs: (hot: ${hot.length})");
    DumpGridStats(grid, opt, evalsBefore, info, form.evals);
  }
  return info;
}

void DoGrid(int pattern, int repeats, double gridPointDistance) {
  print("create grid");
  final dimension = Shapes.kScale / 2.0;
  final GridDimension gd = new GridDimension(
      new VM.Vector3(-dimension, -dimension, -dimension),
      new VM.Vector3(dimension, dimension, dimension),
      gridPointDistance);

  Grid grid = new Grid(gd);
  print("Grid ${grid}");

  Shapes.Form form = allShapes[pattern];
  print("\nFORM: ${form.name}  repeats: ${repeats}");

  DfsOption opt = new DfsOption()
    ..eval = form.evaluator
    ..errorTolerance = 0.0
    ..seqNo = 666;

  print("rounds: ${repeats}");
  Math.Random rng = new Math.Random(1);
  Shapes.Coeffs coeffs = new Shapes.Coeffs(rng, 25, 10.0);

  for (int i = 0; i < repeats; i++) {
    // we animate a little to exercise more code.
    form.Animate(coeffs, i + 0.0);
    MeshInfo golden;
    {
      DateTime begin = new DateTime.now();
      opt.seqNo++;
      golden = AlgoBruteForce(grid, opt, form, true);
      grid.ClearValues();
      DateTime end = new DateTime.now();
      int msec = end.difference(begin).inMilliseconds;
      print("RunningTime: ${msec}ms\n");
    }
    // We get stack overflows otherwise
    if (gridPointDistance >= 20.0) {
      DateTime begin = new DateTime.now();
      opt.seqNo++;
      MeshInfo mi = AlgoRecursive(grid, opt, form, true);
      grid.ClearValues();
      DateTime end = new DateTime.now();
      int msec = end.difference(begin).inMilliseconds;
      print("RunningTime ${msec}ms\n");
      golden.EqualOrDie(mi);
    }

    {
      DateTime begin = new DateTime.now();
      opt.seqNo++;
      MeshInfo mi = AlgoDfs(grid, opt, form, true);
      grid.ClearValues();
      DateTime end = new DateTime.now();
      int msec = end.difference(begin).inMilliseconds;
      print("RunningTime ${msec}ms\n");
      golden.EqualOrDie(mi);
    }
  }
}

int numWorkers = 0;

void GridWorker(List<Object> args) {
  int no = args[0];
  SendPort sendPort = args[1];
  ReceivePort receivePort = new ReceivePort();
  sendPort.send([no, receivePort.sendPort]);
  receivePort.listen((List msg) {
    assert(msg[0] == no);
    print('Received [$msg]');
    DoGrid(msg[1], msg[2], 20.0);
    sendPort.send([no, 'done']);
  });
}

void main(List<String> args) {
  if (args.length == 0) {
    print("invocation: grid_test.dart <form> <repeats> <dist>");
    for (int i = 0; i < allShapes.length; i++) {
      print("${i}:  ${allShapes[i].name}");
    }
    return;
  }
  final int pattern = int.parse(args[0]);
  final int repeats = int.parse(args[1]);
  final double gridPointDistance = double.parse(args[2]);

  if (numWorkers == 0) {
    if (pattern == -1) {
      for (int i = 0; i < allShapes.length; i++) {
        DoGrid(i, repeats, gridPointDistance);
      }
    } else {
      DoGrid(pattern, repeats, gridPointDistance);
    }
    return;
  }

  List<SendPort> sendPort = [];
  for (int i = 0; i < numWorkers; i++) {
    sendPort.add(null);
  }
  // Share for all workers
  ReceivePort receive = new ReceivePort();
  receive.listen((List msg) {
    int no = msg[0];
    if (sendPort[no] == null) {
      sendPort[no] = msg[1];
    } else {
      print('From isolate: $msg');
    }
    sendPort[no].send([no, 1, 1]);
  });

  for (int i = 0; i < numWorkers; i++) {
    Isolate.spawn(GridWorker, [i, receive.sendPort]);
  }
}
