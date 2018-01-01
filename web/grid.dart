library grid;

import 'dart:typed_data';

import 'package:vector_math/vector_math.dart' as VM;
import 'elengyeltables.dart' as tables;

final bool _logTiminig = false;

// LR left-right
// BT bottom-top
// FN far-near

const int LBF = 0; // [x:0,y:0,z:0]
const int LBN = 1;
const int LTF = 2;
const int LTN = 3;
const int RBF = 4;
const int RBN = 5;
const int RTF = 6;
const int RTN = 7; // [x:1,y:1,z:1]

bool threePointNormal(VM.Vector3 a, VM.Vector3 b, VM.Vector3 c, VM.Vector3 tmp,
    VM.Vector3 normal) {
  tmp
    ..setFrom(b)
    ..sub(a);
  normal
    ..setFrom(c)
    ..sub(a);

  normal.crossInto(tmp, normal);
  double len = normal.length;
  if (len == 0) {
    return false;
  }

  return true;
}

void ComputeNormalsFromFaces(
    List<VM.Vector3> vertices, List<VM.Vector3> normals, List<int> faces) {
  assert(normals.length == 0);
  final VM.Vector3 tmp = new VM.Vector3.zero();
  final VM.Vector3 norm = new VM.Vector3.zero();
  for (int i = 0; i < vertices.length; i++) {
    normals.add(new VM.Vector3.zero());
  }
  final int n = faces.length & 0xffffff;
  for (int i = 0; i < n; i += 3) {
    int ia = faces[i + 0];
    int ib = faces[i + 1];
    int ic = faces[i + 2];
    VM.Vector3 a = vertices[ia];
    VM.Vector3 b = vertices[ib];
    VM.Vector3 c = vertices[ic];
    threePointNormal(a, b, c, tmp, norm);
    normals[ia].add(norm);
    normals[ib].add(norm);
    normals[ic].add(norm);
  }
  for (VM.Vector3 v in normals) {
    v.normalize();
  }
}

// This is also available in the shader
void ComputePreciseNormals(
    MeshInfo mesh, SignedDistanceFunction eval, double epsilon) {
  VM.Vector3 v = new VM.Vector3.zero();
  assert(mesh.nNormals3 == 0);
  final int n = mesh.nVertices3 & 0xffffff;
  for (int i = 0; i < n; i += 3) {
    v.x = mesh.vertices[i + 0];
    v.y = mesh.vertices[i + 1];
    v.z = mesh.vertices[i + 2];
    final double val = eval(v);
    v.x -= epsilon;
    double x = eval(v);
    v.x += epsilon;
    v.y -= epsilon;
    double y = eval(v);
    v.y += epsilon;
    v.z -= epsilon;
    double z = eval(v);
    v.z += epsilon;
    v.x = x - val;
    v.y = y - val;
    v.z = z - val;
    v.normalize();
    mesh.AddNormal(v);
  }
  assert(mesh.nNormals3 == mesh.nVertices3);
}

final VM.Vector3 _kZerov3 = new VM.Vector3.zero();

// This represents the surface where the sign flipped and can be
// trivially converted to a MeshData
class MeshInfo {
  Uint32List faces = new Uint32List(3 * 300 * 1024);
  Float32List vertices = new Float32List(3 * 150 * 1024);
  Float32List vertices2 = new Float32List(3 * 150 * 1024);
  Float32List normals = new Float32List(3 * 150 * 1024);
  int nFaces3 = 0;
  int nVertices3 = 0;
  int nNormals3 = 0;

  void AddFace(int a, int b, int c) {
    if (nFaces3 >= faces.length) {
      print("face realloc");
      Uint32List old = faces;
      faces = new Uint32List(old.length * 2);
      for (int i = 0; i < old.length; i++) faces[i] = old[i];
    }
    faces[nFaces3 + 0] = a;
    faces[nFaces3 + 1] = b;
    faces[nFaces3 + 2] = c;
    nFaces3 += 3;
  }

  void AddVertex2(VM.Vector3 v1, VM.Vector3 v2) {
    if (nVertices3 >= vertices.length) {
      print("vertex realloc");
      {
        Float32List old = vertices;
        vertices = new Float32List(old.length * 2);
        for (int i = 0; i < old.length; i++) vertices[i] = old[i];
      }
      {
        Float32List old = vertices2;
        vertices2 = new Float32List(old.length * 2);
        for (int i = 0; i < old.length; i++) vertices2[i] = old[i];
      }
    }
    vertices[nVertices3 + 0] = v1.storage[0];
    vertices[nVertices3 + 1] = v1.storage[1];
    vertices[nVertices3 + 2] = v1.storage[2];

    vertices2[nVertices3 + 0] = v2.storage[0];
    vertices2[nVertices3 + 1] = v2.storage[1];
    vertices2[nVertices3 + 2] = v2.storage[2];
    nVertices3 += 3;
  }

  void AddVertex(VM.Vector3 v) {
    AddVertex2(v, _kZerov3);
  }

  void AddNormal(VM.Vector3 v) {
    if (nNormals3 >= normals.length) {
      print("vertex realloc");
      Float32List old = normals;
      normals = new Float32List(old.length * 2);
      for (int i = 0; i < old.length; i++) normals[i] = old[i];
    }
    normals[nNormals3 + 0] = v.storage[0];
    normals[nNormals3 + 1] = v.storage[1];
    normals[nNormals3 + 2] = v.storage[2];
    nNormals3 += 3;
  }

  Float32List GetVertices() {
    return new Float32List.view(vertices.buffer, 0, nVertices3);
  }

  Float32List GetVertices2() {
    return new Float32List.view(vertices2.buffer, 0, nVertices3);
  }

  Float32List GetNormals() {
    assert(nVertices3 == nNormals3);
    return new Float32List.view(normals.buffer, 0, nNormals3);
  }

  Int32List GetFaces() {
    return new Int32List.view(faces.buffer, 0, nFaces3);
  }

  // For coarse Testing
  void EqualOrDie(MeshInfo o) {
    assert(nFaces3 == o.nFaces3);
    assert(nVertices3 == o.nVertices3);
    assert(nNormals3 == o.nNormals3);
  }

  @override
  String toString() {
    return "MeshInfo: F:${nFaces3/3}  V:${nVertices3/3}";
  }
}

typedef double SignedDistanceFunction(VM.Vector3 v);

class DfsOption {
  // The SDF evaluates distance at given point
  SignedDistanceFunction eval;
  // Euclean neighborhood to explore if we found a sign flip in the distance.
  // This reflects a "continuity" assumption about the SDF which implies
  // that sign flips are more or less adjacent.
  int ttl = 3;
  // Each pass has a unique sequence so that we can ignore results
  // from previous passes.
  int seqNo;

  double errorTolerance = 0.05;
  double epsilon = 0.05;
}

const int POINT_GRID = 1;
const int POINT_CUBE = 2;
const int POINT_LT_ZERO = 4;
const int POINT_GE_ZERO = 8;
const int POINT_ENQUEUED = 16;

int mydiv(double start, double end, double d) {
  assert(start < end);
  int n = 0;
  while (start <= end) {
    start += d;
    n++;
  }
  return n;
}

class GridDimension {
  double d;
  int nx;
  int ny;
  int nz;
  int n;
  VM.Vector3 lbf;
  VM.Vector3 utn;

  GridDimension(this.lbf, VM.Vector3 utnTentative, this.d) {
    nx = mydiv(lbf.x, utnTentative.x, d);
    ny = mydiv(lbf.y, utnTentative.y, d);
    nz = mydiv(lbf.z, utnTentative.z, d);
    utn = new VM.Vector3(
        lbf.x + nx * d - d, lbf.y + ny * d - d, lbf.z + nz * d - d);
    // borders
    nx += 2;
    ny += 2;
    nz += 2;
    n = nx * ny * nz;
  }
}

// This is the workhorse for the cube marching.
class Grid {
  // For each point in the  grid
  // The x,y,z postion
  final Float32List _point_pos;
  // Maybe the  SDF value
  final Float32List _point_value;
  // Maybe 4 ints of reused info
  final Int32List _point_reuse;
  // Maybe the sequence number from DfsOption
  // If this seqNo matches then _point_value has the correct value
  final Uint16List _point_visited_seq;
  // Maybe the sequence number from DfsOption
  final Uint16List _point_completed_seq;

  // Maybe the
  final Uint8List _point_max_ttl;
  // Subset of POINT_BITS from above
  final Uint8List _point_bits;

  //points to process
  final Uint32List _dfs_stack;

  final Uint32List _cornerOffsets = new Uint32List(8);

  // Constants
  final VM.Vector3 _lbf;
  final VM.Vector3 _utn;
  final double _d;
  final int _nx;
  final int _ny;
  final int _nz;
  final int _nextX;
  final int _nextY;
  final int _nextZ;

  // reduce allocations
  final VM.Vector3 _tmp = new VM.Vector3.zero();
  final VM.Vector3 _tmpa = new VM.Vector3.zero();
  final VM.Vector3 _tmpb = new VM.Vector3.zero();

  Grid(GridDimension gd)
      : _lbf = gd.lbf,
        _utn = gd.utn,
        _d = gd.d,
        _nx = gd.nx,
        _ny = gd.ny,
        _nz = gd.nz,
        //
        _nextX = gd.ny * gd.nz,
        _nextY = gd.nz,
        _nextZ = 1,
        //
        _point_reuse = new Int32List(gd.n * 4),
        _point_pos = new Float32List(gd.n * 3),
        _point_value = new Float32List(gd.n),
        _point_visited_seq = new Uint16List(gd.n),
        _point_completed_seq = new Uint16List(gd.n),
        _point_max_ttl = new Uint8List(gd.n),
        _point_bits = new Uint8List(gd.n),
        _dfs_stack = new Uint32List(gd.n) {
    InitArrays();
    //
    _cornerOffsets[0] = 0;
    _cornerOffsets[1] = _nextX;
    _cornerOffsets[2] = _nextY;
    _cornerOffsets[3] = _nextX + _nextY;
    _cornerOffsets[4] = _nextZ;
    _cornerOffsets[5] = _nextZ + _nextX;
    _cornerOffsets[6] = _nextZ + _nextY;
    _cornerOffsets[7] = _nextZ + _nextY + _nextX;
  }

  int GetReusePoint(int c, int direction) {
    if (direction & 1 == 1) c -= _nextX;
    if (direction & 2 == 2) c -= _nextY;
    if (direction & 4 == 4) c -= _nextZ;
    return c;
  }

  @override
  String toString() {
    return "points: ${_nx * _ny * _nz} ${_nx}x${_ny}x${_nz} extremes: ${_lbf} ${_utn}";
  }

  int GetNeighbor(int i, int dx, int dy, int dz) {
    return i + ((dx * _ny + dy) * _nz) + dz;
  }

  bool IsCubeMarchingPoint(int x, int y, int z) {
    return x > 0 && y > 0 && z > 0 && x < _nx - 2 && y < _ny - 2 && z < _nz - 2;
  }

  bool IsRealGridPoint(int x, int y, int z) {
    return x > 0 && y > 0 && z > 0 && x < _nx - 1 && y < _ny - 1 && z < _nz - 1;
  }

  void InitArrays() {
    int i = 0;
    for (int x = 0; x < _nx; x++) {
      for (int y = 0; y < _ny; y++) {
        for (int z = 0; z < _nz; z++) {
          int bits = 0;
          if (IsRealGridPoint(x, y, z)) {
            bits |= POINT_GRID;
            _point_pos[3 * i + 0] = _lbf.x + (x - 1) * _d;
            _point_pos[3 * i + 1] = _lbf.y + (y - 1) * _d;
            _point_pos[3 * i + 2] = _lbf.z + (z - 1) * _d;
          }
          if (IsCubeMarchingPoint(x, y, z)) bits |= POINT_CUBE;
          _point_bits[i] = bits;
          i++;
        }
      }
    }
  }

  int SnapToPoint(VM.Vector3 p) {
    VM.Vector3 o = (p - _lbf) / _d;
    double x = 1.0 + o.x.roundToDouble().clamp(0.0, _nx - 2.0);
    double y = 1.0 + o.y.roundToDouble().clamp(0.0, _ny - 2.0);
    double z = 1.0 + o.z.roundToDouble().clamp(0.0, _nz - 2.0);

    return (x.ceil() * _ny + y.ceil()) * _nz + z.ceil();
  }

  void ClearValues() {
    for (int i = 0; i < _point_value.length; i++) _point_value[i] = 0.0;
  }

  // A point is part of the implicit surface if has pos (neg) distance
  // and ope of its six neighbors has neg (pos) distance.
  List<int> GetSurfacePoints(bool below, DfsOption opt) {
    List<int> out = [];
    final int n = _nx * _ny * _nz;
    bool IsMatch(int x, bool below) {
      if (_point_visited_seq[x] != opt.seqNo) return false;
      if (_point_bits[x] & POINT_GRID == 0) return false;
      return below == (_point_value[x] < 0.0);
    }

    for (int i = 0; i < n; i++) {
      if (!IsMatch(i, below)) continue;
      if (IsMatch(i + _nextX, !below) ||
          IsMatch(i - _nextX, !below) ||
          IsMatch(i + _nextY, !below) ||
          IsMatch(i - _nextY, !below) ||
          IsMatch(i + _nextZ, !below) ||
          IsMatch(i - _nextZ, !below)) {
        out.add(i);
      }
    }

    return out;
  }

  double UpdateValue(int p, SignedDistanceFunction eval, int seqNo) {
    _tmp.copyFromArray(_point_pos, p * 3);
    final double val = eval(_tmp);
    _point_value[p] = val;
    _point_visited_seq[p] = seqNo;
    _point_completed_seq[p] = seqNo;
    _point_reuse[4 * p + 0] = -1;
    _point_reuse[4 * p + 1] = -1;
    _point_reuse[4 * p + 2] = -1;
    _point_reuse[4 * p + 3] = -1;
    return val;
  }

  // Recompute values for all grid points
  void ValueUpdateBruteForce(final DfsOption opt) {
    final int seqNo = opt.seqNo;
    final SignedDistanceFunction eval = opt.eval;
    DateTime begin = new DateTime.now();
    final int n = _nx * _ny * _nz;
    for (int i = 0; i < n; i++) {
      int currentBits = _point_bits[i];
      if (currentBits & POINT_GRID == 0) continue;
      _point_max_ttl[i] = opt.ttl;
      double value = UpdateValue(i, eval, seqNo);
      if (value < 0.0) {
        currentBits |= POINT_LT_ZERO;
      } else {
        currentBits &= ~POINT_LT_ZERO;
      }
      _point_bits[i] = currentBits;
    }
    DateTime end = new DateTime.now();
    int msec = end.difference(begin).inMilliseconds;
    if (_logTiminig) print("ValueUpdateBruteForce [${msec}ms]");
  }

  // TODO: move this into the vertex shader
  void InterpolateCorners(int p1, int p2, VM.Vector3 out,
      SignedDistanceFunction eval, double tolerance) {
    double v1 = _point_value[p1];
    double v2 = _point_value[p2];
    _tmpa.copyFromArray(_point_pos, 3 * p1);
    _tmpb.copyFromArray(_point_pos, 3 * p2);

    for (int i = 0; i < 3; ++i) {
      assert(v1 < 0.0);
      assert(v2 >= 0);
      final d = v1 / (v1 - v2);
      assert(1 >= d && d >= 0);
      VM.Vector3.mix(_tmpa, _tmpb, d, out);
      if (tolerance <= 0.0) return;
      if ((v1 - v2).abs() < 2.0 * tolerance) return;
      final double v = eval(out);
      if (v.abs() < tolerance) return;
      if (v < 0.0) {
        v1 = v;
        _tmpa.setFrom(out);
      } else {
        v2 = v;
        _tmpb.setFrom(out);
      }
    }
  }

  int MarchOneCubeInit(int c, int seqNo) {
    int mask = 0;
    // TODO: shouldn't we confirm that the values is valid?
    for (int i = 0; i < 8; i++) {
      final int p = c + _cornerOffsets[i];
      // avoid "incomplete" cubes
      if (seqNo != _point_visited_seq[p]) {
        // This seems to be ok - explain
        // print ("@@@@@ missing value");
        return 0;
      }
      if ((_point_bits[p] & POINT_LT_ZERO) != 0) {
        // the and'ing seems to prevent boxing and has a measurable impact
        mask |= (1 << i) & 0xff;
      }
    }
    // if there are no differences there is nothing to do
    if (mask == 0 || mask == 255) {
      //print ("@@@@@@@@ ${mask}");
      return 0;
    }
    return mask;
  }

  void MarchOneCubeComputeVertices(int c, int seqNo, int mask,
      Uint32List cellVerticesIndices, MeshInfo mesh, final DfsOption opt) {
    final double tolerance = opt.errorTolerance;
    final SignedDistanceFunction eval = opt.eval;
    final int start = mask * tables.vertexDataRowLength;
    for (int i = 0; i < tables.vertexDataRowLength; i++) {
      final int v = tables.vertexData[start + i];
      if (v == 0) break;
      // 1, 2 or 3
      final int reuseSlot = (v >> 8) & 0xf;
      assert(reuseSlot < 4);
      final int reuseDir = (v >> 12) & 0xf;
      int reusePoint = c;
      if (reuseDir != 8) {
        int p2 = GetReusePoint(c, reuseDir);
        if (_point_visited_seq[p2] != seqNo) {
          print("Rare point re-compute required");
          _point_visited_seq[p2] = seqNo;
          _point_reuse[4 * p2 + 0] = -1;
          _point_reuse[4 * p2 + 1] = -1;
          _point_reuse[4 * p2 + 2] = -1;
          _point_reuse[4 * p2 + 3] = -1;
        }
        reusePoint = p2;
      }

      int index = _point_reuse[4 * reusePoint + reuseSlot];
      if (index == -1) {
        final int p1 = c + _cornerOffsets[v & 0xf];
        final int p2 = c + _cornerOffsets[(v >> 4) & 0xf];
        index = mesh.nVertices3 ~/ 3;
        final bool swap = mask & (1 << (v & 0xf)) == 0;

        if (swap) {
          InterpolateCorners(p2, p1, _tmp, eval, tolerance);
        } else {
          InterpolateCorners(p1, p2, _tmp, eval, tolerance);
        }
        mesh.AddVertex(_tmp);
        /*
        _tmpa.copyFromArray(_point_pos, 3 * p1);
        _tmpb.copyFromArray(_point_pos, 3 * p2);
        if (swap) {
          mesh.AddVertex2(_tmpb, _tmpa);
        } else {
          mesh.AddVertex2(_tmpa, _tmpb);
        }
        */
        _point_reuse[4 * reusePoint + reuseSlot] = index;
      }
      cellVerticesIndices[i] = index;
    }
  }

  // If UpdatePoint return true, extract the triangles vertices into v.
  // and the corresponding normals into n

  void MarchOneCubeComputeFaces(
      int mask, Uint32List cellVerticesIndices, MeshInfo mesh) {
    final int start = tables.hotCornerMap[mask] * tables.configClassRowLength;
    final int end = start + 2 + tables.configClass[start + 1] * 3;
    for (int i = start + 2; i < end; i += 3) {
      int ia = cellVerticesIndices[tables.configClass[i + 0]];
      int ib = cellVerticesIndices[tables.configClass[i + 1]];
      int ic = cellVerticesIndices[tables.configClass[i + 2]];
      mesh.AddFace(ia, ib, ic);
    }
  }

  // Each point in hot is consider to be a cube to be processed
  MeshInfo CubeMarch(
      Iterable<int> hot, final DfsOption opt, bool computeNormals) {
    final SignedDistanceFunction eval = opt.eval;
    final double epsilon = opt.epsilon;
    final int seqNo = opt.seqNo;
    MeshInfo mesh = new MeshInfo();
    DateTime begin = new DateTime.now();
    int action = 0;
    int visited = 0;
    // The eight corner points of the cube
    // face vertices for the current cube - at most 12
    final Uint32List cellVerticesIndices = new Uint32List(12);
    for (int i in hot) {
      // assert(_point_bits[i] & POINT_CUBE != 0);
      // assert(_point_visited_seq[i] == seqNo); // we used to ignore these
      visited++;
      final int mask = MarchOneCubeInit(i, seqNo);
      if (mask == 0) continue;
      action++;
      MarchOneCubeComputeVertices(
          i, seqNo, mask, cellVerticesIndices, mesh, opt);
      MarchOneCubeComputeFaces(mask, cellVerticesIndices, mesh);
    }

    if (computeNormals) {
      ComputePreciseNormals(mesh, eval, epsilon);
    }

    DateTime end = new DateTime.now();
    int msec = end.difference(begin).inMilliseconds;
    if (_logTiminig)
      print(
          "CubeMarch [${msec}ms]: visited: ${visited} hotcubes: ${action} vertices: ${mesh.nVertices3 ~/ 3} faces: ${mesh.nFaces3 ~/ 3}");
    return mesh;
  }

  // Give a location inside the object find a surface point
  int FindSurfacePointFromSeedBits(
      VM.Vector3 start, SignedDistanceFunction eval) {
    final int p = SnapToPoint(start);

    //  assert(_point_bits[p] & POINT_GRID != 0);
    if (_point_bits[p] & POINT_GRID == 0) {
      print("interior point not really inside ${start}");
      return -1;
    }
    _tmp.copyFromArray(_point_pos, 3 * p);
    double value = eval(_tmp);
    final bool initial = value < 0;

    for (int dir in [-_nextX, _nextX, -_nextY, _nextY, -_nextZ, _nextZ]) {
      int q = p;
      while (true) {
        q += dir;
        if (_point_bits[q] & POINT_GRID == 0) break;
        _tmp.copyFromArray(_point_pos, 3 * q);
        double value = eval(_tmp);
        if (initial != (value < 0)) return q;
      }
    }

    //assert(false);
    return -1;
  }

  void ValueUpdateRecursive(Iterable<int> points, int ttl, bool wantBelow,
      List<int> hot, final DfsOption opt) {
    final int seqNo = opt.seqNo;
    final SignedDistanceFunction eval = opt.eval;
    DateTime begin = new DateTime.now();
    //print ("${i} ${ttl} ${below}");
    bool dfs(int p, int ttl, bool wantBelow) {
      assert(ttl >= 0);
      int currentBits = _point_bits[p];
      if (currentBits & POINT_GRID == 0) return false;

      if (_point_visited_seq[p] != seqNo) {
        double value = UpdateValue(p, eval, seqNo);
        if (value < 0.0) {
          currentBits |= POINT_LT_ZERO;
        } else {
          currentBits &= ~POINT_LT_ZERO;
        }
        _point_bits[p] = currentBits;
        _point_max_ttl[p] = 0;
        hot.add(p);
      }

      final bool isBelow = _point_value[p] < 0.0;
      if (wantBelow == isBelow) {
        ttl = opt.ttl;
      }
      // we have explored this point prior with a higher ttl.
      if (ttl <= _point_max_ttl[p]) return wantBelow == isBelow;
      _point_max_ttl[p] = ttl;

      if (ttl <= 0) return isBelow == wantBelow;
      final int t = ttl - 1;
      dfs(p - _nextX, t, !isBelow);
      dfs(p + _nextX, t, !isBelow);
      dfs(p - _nextY, t, !isBelow);
      dfs(p + _nextY, t, !isBelow);
      dfs(p - _nextZ, t, !isBelow);
      dfs(p + _nextZ, t, !isBelow);
      return isBelow == wantBelow;
    }

    for (int p in points) {
      dfs(p, ttl, wantBelow);
    }
    DateTime end = new DateTime.now();
    int msec = end.difference(begin).inMilliseconds;
    if (_logTiminig) print("ValueUpdateRecursive [${msec}ms] ${hot.length}");
  }

  // Standard DFS using a queue starting from one or more seedpoints
  void ValueUpdateDfs(Iterable<int> points, bool xwantBelow, List<int> hot,
      final DfsOption opt) {
    final int seqNo = opt.seqNo;
    final SignedDistanceFunction eval = opt.eval;
    final int xttl = opt.ttl;
    DateTime begin = new DateTime.now();
    int sp = 0;

    //int count = 0;
    void enqueue(int p, int ttl, bool wantBelow) {
      int currentBits = _point_bits[p];
      // ignore points at the border
      if (currentBits & POINT_GRID == 0) return;

      int maxTtl;
      if (_point_visited_seq[p] != seqNo) {
        double value = UpdateValue(p, eval, seqNo);
        if (value < 0.0) {
          currentBits |= POINT_LT_ZERO;
        } else {
          currentBits &= ~POINT_LT_ZERO;
        }
        _point_bits[p] = currentBits;
        maxTtl = 0;
      } else {
        maxTtl = _point_max_ttl[p];
      }

      final bool isBelow = (currentBits & POINT_LT_ZERO) != 0;
      if (wantBelow == isBelow) {
        // we found a transition - start again with a full ttl budget
        ttl = xttl;
      }
      // we have seen this point before with more range
      if (ttl <= maxTtl) return;
      // This is the first time we process this transition point
      if (maxTtl == 0 && ttl > 0 && (currentBits & POINT_CUBE) != 0) hot.add(p);
      _point_max_ttl[p] = ttl;
      if (ttl <= 1) return;

      if ((currentBits & POINT_ENQUEUED) != 0) return;
      // count++;
      _point_bits[p] = currentBits | POINT_ENQUEUED;
      _dfs_stack[sp] = p;
      sp++;
    }

    for (int p in points) {
      enqueue(p, xttl, xwantBelow);
    }

    // invariants: everything in the stack has been evaluated
    while (sp > 0) {
      //print ("${stack.length}");
      sp--;
      final int p = _dfs_stack[sp];
      int bits = _point_bits[p];
      _point_bits[p] = bits & ~POINT_ENQUEUED;
      final int t = _point_max_ttl[p] - 1;
      final bool isBelow = (bits & POINT_LT_ZERO) != 0;
      // stuff with ttl = 0 gets evaluated but not enqueued
      enqueue(p - _nextX, t, !isBelow);
      enqueue(p + _nextX, t, !isBelow);
      enqueue(p - _nextY, t, !isBelow);
      enqueue(p + _nextY, t, !isBelow);
      enqueue(p - _nextZ, t, !isBelow);
      enqueue(p + _nextZ, t, !isBelow);
    }

    DateTime end = new DateTime.now();
    int msec = end.difference(begin).inMilliseconds;
    if (_logTiminig) print("ValueUpdateDfs [${msec}ms] ${hot.length}");
  }

  List<int> ValidCubePoints() {
    List<int> out = <int>[];
    final int n = _nx * _ny * _nz;
    for (int i = 0; i < n; i++) {
      if ((_point_bits[i] & POINT_CUBE) != 0) out.add(i);
    }

    return out;
  }
}
