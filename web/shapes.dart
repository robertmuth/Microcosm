library shapes;

import 'dart:core';
import 'dart:math' as Math;
import 'dart:typed_data';

import 'package:vector_math/vector_math.dart' as VM;

const double epsilon = 0.000000001;
const double kInfinity = 10000.0;
const double TWOPI = 2.0 * Math.PI;
const double HALFPI = 0.5 * Math.PI;
const double kScale = 1000.0;

final VM.Vector3 zero3 = new VM.Vector3.zero();
final VM.Vector3 e011 = new VM.Vector3(0.0, 1.0, 1.0);

double length2(double x, double y, double z) {
  return x * x + y * y + z * z;
}

// ============================================================
// ============================================================
abstract class IsoShape {
  // TODO: make this final
  VM.Matrix4 transform = new VM.Matrix4.identity();
  final VM.Matrix4 inv_transform = new VM.Matrix4.identity();
  final VM.Vector4 params = new VM.Vector4.zero();
  // use for eval
  final VM.Vector3 vtmp = new VM.Vector3.zero();

  // Returns > 1 out side and < 1 inside.
  // Must not modify pos
  // Note, this is unusual and you may want to use the inverse instead.
  //
  double EvalAt(VM.Vector3 pos);

  String CodeForEval() {
    assert(false);
    return "";
  }

  void ComputeInvTransform() {
    inv_transform.copyInverse(transform);
  }

  VM.Vector3 InteriorPoint() {
    return transform * kOrigin;
  }
}

final VM.Vector3 kOrigin = new VM.Vector3.zero();

class Sphere extends IsoShape {
  final double _radius2_inv;

  Sphere(double radius) : _radius2_inv = 1.0 / (radius * radius) {
    params[0] = _radius2_inv;
  }

  @override
  double EvalAt(VM.Vector3 pos) {
    inv_transform.transformed3(pos, vtmp);
    return vtmp.length2 * _radius2_inv;
  }

  @override
  String CodeForEval() {
    return "SDF1Sphere(p, param.x);";
  }
}

class SphereUnit extends IsoShape {
  SphereUnit();

  @override
  double EvalAt(VM.Vector3 pos) {
    inv_transform.transformed3(pos, vtmp);
    return vtmp.length2;
  }

  @override
  String CodeForEval() {
    return "SDF1SphereUnit(p);";
  }
}

class Capsule extends IsoShape {
  double _radius2_inv;
  double _halflength;
  double _radius_orig;
  double _halflength_orig;

  Capsule(this._radius_orig, this._halflength_orig)
      : _radius2_inv = 1.0 / (_radius_orig * _radius_orig),
        _halflength = _halflength_orig {
    params[0] = _radius2_inv;
    params[1] = _halflength;
  }

  void SetScale2(double sxy, double sz) {
    if (sxy == 0.0 || sz == 0) throw "scale may not be zero";
    _radius2_inv = 1.0 / (_radius_orig * _radius_orig * sxy * sxy);
    _halflength = _halflength_orig * sz;
    params[0] = _radius2_inv;
    params[1] = _halflength;
  }

  @override
  double EvalAt(VM.Vector3 pos) {
    inv_transform.transformed3(pos, vtmp);

    double z = vtmp.z.abs() - _halflength;
    return length2(vtmp.x, vtmp.y, z > 0.0 ? z : 0.0) * _radius2_inv;
  }

  @override
  String CodeForEval() {
    return "SDF1Capsule(p, param.x, param.y);";
  }
}

class Torus extends IsoShape {
  final VM.Vector3 _inner = new VM.Vector3.zero();
  double _thickness2_inv;
  double _orbit;
  double _thickness_orig;
  double _orbit_orig;

  // orbit is the center radius of the ring
  // thickness is radius of sausage wrapped in a circle
  Torus(this._thickness_orig, this._orbit_orig)
      : _thickness2_inv = 1.0 / (_thickness_orig * _thickness_orig),
        _orbit = _orbit_orig {
    _inner.x = _orbit;
    params[0] = _orbit;
    params[1] = _thickness2_inv;
  }

  void SetScale2(double st, double so) {
    if (st == 0.0 || so == 0.0) throw "scale may not be zero";
    _thickness2_inv = 1.0 / (_thickness_orig * _thickness_orig * st * st);
    _orbit = _orbit_orig * so;
    _inner.x = _orbit;
    params[0] = _orbit;
    params[1] = _thickness2_inv;
  }

  @override
  double EvalAt(VM.Vector3 pos) {
    inv_transform.transformed3(pos, vtmp);
    double t = Math.sqrt(vtmp.x * vtmp.x + vtmp.y * vtmp.y) - _orbit;
    return (t * t + vtmp.z * vtmp.z) * _thickness2_inv;
  }

  @override
  String CodeForEval() {
    return "SDF1Torus(p, param.x, param.y);";
  }

  @override
  VM.Vector3 InteriorPoint() {
    return transform * _inner;
  }
}

class Cube extends IsoShape {
  double _inv_rx2;
  double _inv_ry2;
  double _inv_rz2;

  Cube(double rx, double ry, double rz)
      : _inv_rx2 = 1 / (rx * rx),
        _inv_ry2 = 1 / (ry * ry),
        _inv_rz2 = 1 / (rz * rz) {
    params[0] = _inv_rx2;
    params[1] = _inv_ry2;
    params[2] = _inv_rz2;
  }

  // Returns > 1 out side and < 1 inside.
  @override
  double EvalAt(VM.Vector3 pos) {
    inv_transform.transformed3(pos, vtmp);
    double dx = (vtmp.x * vtmp.x) * _inv_rx2;
    double dy = (vtmp.y * vtmp.y) * _inv_ry2;
    double dz = (vtmp.z * vtmp.z) * _inv_rz2;
    return Math.max(Math.max(dx, dy), dz);
  }

  @override
  String CodeForEval() {
    return "SDF1Cube(p, param.x, param.y, param.z);";
  }
}

// Cube like
class RoundedCube extends IsoShape {
  final double _rx;
  final double _ry;
  final double _rz;
  final double _e2;

  RoundedCube(double edge, this._rx, this._ry, this._rz)
      : _e2 = 1.0 / (edge * edge) {
    params[0] = _rx;
    params[1] = _ry;
    params[2] = _rz;
    params[3] = _e2;
  }

  // Returns > 1 out side and < 1 inside.
  @override
  double EvalAt(VM.Vector3 pos) {
    inv_transform.transformed3(pos, vtmp);

    double z = vtmp.z.abs() - _rz;
    double y = vtmp.y.abs() - _ry;
    double x = vtmp.x.abs() - _rx;
    double l2 =
        length2(x > 0.0 ? x : 0.0, y > 0.0 ? y : 0.0, z > 0.0 ? z : 0.0);
    return l2 * _e2;
  }

  @override
  String CodeForEval() {
    return "SDF1RoundedCube(p, param.x, param.y, param.z, params.w);";
  }
}

class RoundedCubeUnit extends IsoShape {
  RoundedCubeUnit();

  // Returns > 1 out side and < 1 inside.
  @override
  double EvalAt(VM.Vector3 pos) {
    inv_transform.transformed3(pos, vtmp);

    double z = vtmp.z.abs() - 1.0;
    double y = vtmp.y.abs() - 1.0;
    double x = vtmp.x.abs() - 1.0;
    double l2 =
        length2(x > 0.0 ? x : 0.0, y > 0.0 ? y : 0.0, z > 0.0 ? z : 0.0);
    return l2 * 4.0;
  }

  @override
  String CodeForEval() {
    return "SDF1RoundedCubeUnit(p);";
  }
}

class Knot extends IsoShape {
  final VM.Vector3 _center = new VM.Vector3.zero();
  final double _r1_orig;
  final double _r2_orig;
  final double _thickness2;
  final int _coils;
  final int _twists;
  final double _lat_offset;
  final double _twists_over_coils;
  double _r1;
  double _r2;

  // standard coils = 3, twists = 2, r1=1000, r2=500
  Knot(
      this._r1_orig, this._r2_orig, double thickness, this._coils, this._twists)
      : _lat_offset = TWOPI / _coils,
        _twists_over_coils = _twists / _coils,
        _thickness2 = thickness * thickness,
        _r1 = _r1_orig,
        _r2 = _r2_orig {
    params[0] = _r1;
    params[1] = _r2;
  }

  void SetScale2(double s1, double s2) {
    if (s1 == 0.0 || s2 == 0.0) throw "scale may not be zero";
    _r1 = _r1_orig * s1;
    _r2 = _r2_orig * s2;
    params[0] = _r1;
    params[1] = _r2;
  }

  // Returns > 1 out side and < 1 inside.
  @override
  double EvalAt(VM.Vector3 pos) {
    inv_transform.transformed3(pos, vtmp);
    final double tmp = Math.sqrt(vtmp.x * vtmp.x + vtmp.y * vtmp.y) - _r1;
    final double lat = Math.atan2(vtmp.y, vtmp.x) * _twists_over_coils;
    double r = 0.0;
    for (int i = 0; i < _coils; i++) {
      final double lon = lat + _lat_offset * i;
      final double hor = tmp - Math.cos(lon) * _r2;
      final double ver = vtmp.z - Math.sin(lon) * _r2;
      final double d = hor * hor + ver * ver;
      if (d == 0.0) {
        return 0.0;
      }
      r += _thickness2 / d;
    }
    if (r == 0.0) return kInfinity;
    return 1.0 / r;
  }

  @override
  String CodeForEval() {
    return "SDF1Knot(p, ${_coils}, ${_twists_over_coils}, ${_lat_offset}, ${_thickness2.toStringAsFixed(2)}, param.x, param.y)";
  }

  @override
  VM.Vector3 InteriorPoint() {
    // The original had serveral seed points
    // for(int i=0; i<coils; ++i){
    // _center.x = _r1 + Math.cos( i * _lat_offset) * _r2;
    // _center.z = Math.sin(i * _lat_offset) * _r2;
    // }
    _center.x = _r1 + _r2;
    return _center;
  }
}

class Heart extends IsoShape {
  final double _r;

  Heart(this._r) {
    params[0] = _r;
  }

// Returns > 1 out side and < 1 inside.
// For more accuracy we could apply the sqrt
  @override
  double EvalAt(VM.Vector3 pos) {
    inv_transform.transformed3(pos, vtmp);
    final double x = vtmp.x / _r * 2.0;
    final double y = vtmp.y / _r;
    final double z = vtmp.z / _r;
    final double xx = x * x;
    final double yy = y * y;
    final double zz = z * z;
    final double zzz = zz * z;

    double a = 2 * xx + yy + zz - 1;
    a = a * a * a;
    double b = 0.1 * xx * zzz + yy * zzz;
    return (a - b) + 1;
  }

  @override
  String CodeForEval() {
    return "SDF1Heart(p, param.x)";
  }
}

// ============================================================
// ============================================================

class Coeffs {
  final int num;
  final Float32List _rateSeed;
  final Float32List _phaseSeed;

  Coeffs(Math.Random rng, this.num, double speed)
      : _rateSeed = new Float32List(num),
        _phaseSeed = new Float32List(num) {
    for (int i = 0; i < num; i++) {
      _rateSeed[i] = rng.nextDouble() * 0.015 * speed;
      _phaseSeed[i] = rng.nextDouble() * TWOPI - Math.PI;
    }
  }

  Float32List get rateSeed => _rateSeed;
  Float32List get phaseSeed => _phaseSeed;

  double phase(int i, double now) {
    return (_phaseSeed[i] + _rateSeed[i] * now).remainder(TWOPI);
  }

  double value(int i, double now, [double freq = 1.0]) {
    return Math.cos(phase(i, now) * freq);
  }
}

// ============================================================
// ============================================================

abstract class Form {
  int evals = 0;
  final String name;
  List<IsoShape> _shapes = [];

  Form(this.name);

  // helper
  VM.Matrix4 _r = new VM.Matrix4.zero();

  VM.Matrix4 makeRotX(double angle) {
    return _r
      ..setIdentity()
      ..setRotationX(angle);
  }

  VM.Matrix4 makeRotY(double angle) {
    return _r
      ..setIdentity()
      ..setRotationY(angle);
  }

  VM.Matrix4 makeRotZ(double angle) {
    return _r
      ..setIdentity()
      ..setRotationZ(angle);
  }

  VM.Matrix4 makeTrans(double x, double y, double z) {
    return _r
      ..setIdentity()
      ..setTranslationRaw(x, y, z);
  }

  VM.Matrix4 make9and12(double s, double t) {
    _r.setIdentity();
    _r[9] = s;
    _r[12] = t;
    return _r;
  }

  VM.Matrix4 makeScale(double x, double y, double z) {
    _r.setIdentity();
    _r[0] = x;
    _r[5] = y;
    _r[10] = z;
    return _r;
  }

  void setTrans(VM.Matrix4 m, double x, double y, double z) {
    m
      ..setIdentity()
      ..setTranslationRaw(x, y, z);
  }

  void addShape(IsoShape s) {
    _shapes.add(s);
  }

  // The individual shapes return
  // returns > 1 outside and < 1 inside.
  // E.g. for a sphere we get d^2 / r^2
  //       (points close to the center have small values)/
  // This function inverts that:
  // < 1 means outside, > 1 means inside.
  // Furthermore we standardize on zero iso surface
  // < 0 means outside, > 0 means inside.
  double evaluator(VM.Vector3 pos) {
    evals++;
    double value = 0.0;
    // for loop is faster than using iterator
    for (int i = 0; i < _shapes.length; ++i) {
      // This cause an allocation:
      // VM.Vector3 new_pos = _inv_transforms[i] * pos;
      // print("${pos} -> ${new_pos}");
      // double v = _shapes[i].EvalAt(new_pos);
      // This tries to do the same without
      double v = _shapes[i].EvalAt(pos);
      //print("eval ${pos}[$i]: ${value}");
      value += v == 0.0 ? kInfinity : 1.0 / v;
    }
    //print("eval ${pos}: ${value}");
    return value - 1.0;
  }

  Iterable<VM.Vector3> InteriorPoints() {
    List<VM.Vector3> out = new List<VM.Vector3>(_shapes.length);
    for (int i = 0; i < _shapes.length; i++) {
      out[i] = _shapes[i].InteriorPoint();
    }
    return out;
  }

  void UpdateShapeInvTransform(Float32List out) {
    assert(16 * _shapes.length < out.length);
    for (int i = 0; i < _shapes.length; i++) {
      final int o = 16 * i;
      Float32List src = _shapes[i].inv_transform.storage;
      for (int j = 0; j < 16; ++j) {
        out[o + j] = src[j];
      }
    }
  }

  void UpdateShapeParams(Float32List out) {
    assert(4 * _shapes.length < out.length);
    for (int i = 0; i < _shapes.length; i++) {
      final int o = 4 * i;
      Float32List src = _shapes[i].params.storage;
      for (int j = 0; j < 4; ++j) {
        out[o + j] = src[j];
      }
    }
  }

  String ExtractShaderCodeForEval() {
    List<String> out = <String>[];
    out.add("""
float FormEval(vec3 pos) {
    float res = 0.0;
""");

    for (int i = 0; i < _shapes.length; i++) {
      out.add("""
    {
        vec3 p = (uShapeInvTransforms[$i] * vec4(pos, 1.0)).xyz;
        vec4 param = uShapeParams[$i];
        float v = ${_shapes[i].CodeForEval()};
        if (v == 0.0) return 10000.0;
        res += 1.0 / v;
    }
""");
    }
    out.add("""
    return res;
}
""");
    return out.join("");
  }

  // must be overloaded
  void Animate(Coeffs coeffs, double now);
}

class Metaballs extends Form {
  final double sphereScale = 0.06 * kScale;
  List<SphereUnit> spheres = <SphereUnit>[];

  // 5 is a good  number
  Metaballs(int num) : super("Metaballs(${num})") {
    for (int i = 0; i < num; i++) {
      SphereUnit s = new SphereUnit();
      addShape(s);
      spheres.add(s);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    double x = 0.0;
    final double dx = TWOPI / _shapes.length;

    for (SphereUnit s in spheres) {
      VM.Matrix4 m = s.transform;
      m.setIdentity();
      m.rotate(e011, coeffs.phase(10, now) * 3.0 + x);
      m.multiply(makeRotX(coeffs.phase(9, now) * 2.0 + x));
      m.multiply(makeTrans(0.0, 0.27 * kScale * coeffs.value(7, now),
          0.27 * kScale * coeffs.value(8, now)));
//      m.multiply(makeScale(sphereScale, sphereScale, sphereScale));
//      m.setIdentity();
//      // first move  (rotating a sphere is a nop in out setup)
//      m.setTranslationRaw(0.0, 0.27 * kScale * coeffs.value(7, now),
//          0.27 * kScale * coeffs.value(8, now));
//      // cool effect
//      //m = m * xrotate(e100, coeffs.phase(9, now) * 2.0 + x);
//      m = makeRot(e100, coeffs.phase(9, now) * 2.0 + x) * m;
//      m = makeRot(e011, coeffs.phase(10, now) * 3.0 + x) * m;
//
      for (int i = 0; i < 12; i++) {
        m[i] = 0.0;
      }
      m[0] = sphereScale;
      m[5] = sphereScale;
      m[10] = sphereScale;
      m[15] = 1.0;
      s.ComputeInvTransform();
      x += dx;
    }
  }
}

class TriangleOfSpheres extends Form {
  final double sphereScale = 0.06 * kScale;
  final double sphereDist = 0.28 * kScale;
  final List<SphereUnit> _spheres = <SphereUnit>[];

  // 5, 7, 8 are good  numbers
  TriangleOfSpheres(int num) : super("TriangleOfSpheres(${num})") {
    for (int i = 0; i < num; i++) {
      SphereUnit s = new SphereUnit();
      addShape(s);
      _spheres.add(s);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    double x = coeffs.phase(9, now) * 4.0;
    final double dx = TWOPI / _shapes.length;

    for (SphereUnit s in _spheres) {
      VM.Matrix4 m = s.transform;
      m.setIdentity();
      m.setRotationZ(x);
      m.multiply(makeTrans(
          sphereDist * Math.cos(coeffs.phase(8, now) + HALFPI), 0.0, 0.0));
      m.multiply(makeRotZ(x * -3.0));
      m.multiply(
          makeTrans(sphereDist * Math.cos(coeffs.phase(8, now)), 0.0, 0.0));

      // just preserve the translation
      for (int i = 0; i < 12; i++) {
        m[i] = 0.0;
      }
      m[0] = sphereScale;
      m[5] = sphereScale;
      m[10] = sphereScale;
      m[15] = 1.0;
      //m.scale(sphereScale, sphereScale, sphereScale);

      s.ComputeInvTransform();
      x += dx;
    }
  }
}

class StringOfEllipsoids extends Form {
  final double _offset;

  final double sphereScale = 0.04 * kScale;
  final List<SphereUnit> _spheres = <SphereUnit>[];

  StringOfEllipsoids(int num, this._offset)
      : super("StringOfEllipsoids(${num}, ${_offset})") {
    for (int i = 0; i < num; i++) {
      SphereUnit s = new SphereUnit();
      addShape(s);
      _spheres.add(s);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    double n = 0.0;
    for (SphereUnit s in _spheres) {
      double tr(double a, double b) {
        return Math.cos(a * 3.0 + n * (1.0 + 0.5 * b)) * 0.3 * kScale;
      }

      double sc(double ph) {
        return sphereScale * (Math.cos(ph * 4.0 + n) * 0.4 + 1.2);
      }

      VM.Matrix4 m = s.transform;
      m.setIdentity();
      // fun variation: only change one dimension
      m.setTranslationRaw(
          tr(coeffs.phase(13, now), coeffs.value(14, now)),
          tr(coeffs.phase(15, now), coeffs.value(16, now)),
          tr(coeffs.phase(17, now), coeffs.value(18, now)));
      m.multiply(makeRotZ(Math.cos(coeffs.phase(12, now) * 2.0 + n)));
      m.multiply(makeRotY(Math.cos(coeffs.phase(11, now) * 2.0 + n)));
      m.multiply(makeRotX(Math.cos(coeffs.phase(10, now) * 2.0 + n)));
      m.multiply(makeScale(sc(coeffs.phase(7, now)), sc(coeffs.phase(8, now)),
          sc(coeffs.phase(9, now))));

      s.ComputeInvTransform();

      n += _offset;
    }
  }
}

class Brain extends Form {
  final double sphereScale = 0.05 * kScale;
  final List<SphereUnit> _spheres = <SphereUnit>[];

  // Three clusters of metaballs each of size num
  // 3, 4, 5, 6
  Brain(int num) : super("Brain(${num})") {
    for (int i = 0; i < num * 3; i++) {
      SphereUnit s = new SphereUnit();
      addShape(s);
      _spheres.add(s);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    double x = 0.0;
    final double dx = TWOPI / (_shapes.length / 3);

    for (int i = 0; i < _spheres.length; i += 3) {
      final double offset1 = x + coeffs.phase(9, now) * 4.0;
      final double scale = (1.0 - Math.cos(offset1).abs()) * 2.0 + 1.0;

      final double offset2 = x + coeffs.phase(9, now);
      final double cos1 = Math.cos(offset1);
      final double cos2 = Math.cos(offset2);
      final double sin2 = Math.sin(offset2);
      final double mult = 0.38 * kScale;

      {
        VM.Matrix4 m = _spheres[i + 0].transform;
        setTrans(m, cos1 * mult, cos2 * mult, sin2 * mult);
        m.multiply(makeScale(sphereScale * scale, sphereScale, sphereScale));
        _spheres[i + 0].ComputeInvTransform();
      }
      {
        VM.Matrix4 m = _spheres[i + 1].transform;
        setTrans(m, cos2 * mult, cos1 * mult, sin2 * mult);
        m.multiply(makeScale(sphereScale, sphereScale * scale, sphereScale));
        _spheres[i + 1].ComputeInvTransform();
      }
      {
        VM.Matrix4 m = _spheres[i + 2].transform;
        setTrans(m, cos2 * mult, sin2 * mult, cos1 * mult);
        m.multiply(makeScale(sphereScale, sphereScale, sphereScale * scale));
        _spheres[i + 2].ComputeInvTransform();
      }
      x += dx;
    }
  }
}

class Flower extends Form {
  final double _offset;
  final double sphereScale = 0.1 * kScale;

  final List<SphereUnit> _spheres = <SphereUnit>[];
  // 10, 14, 18
  Flower(int num, this._offset) : super("Flower(${num}, ${_offset})") {
    for (int i = 0; i < num; i++) {
      SphereUnit s = new SphereUnit();
      addShape(s);
      _spheres.add(s);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    double x = 0.0;
    double n = 0.0;
    double dx = TWOPI / _shapes.length;

    for (SphereUnit s in _spheres) {
      final double z =
          Math.cos(coeffs.phase(0, now) * 3.0 + n * coeffs.value(1, now));
      VM.Matrix4 m = s.transform;
      m.setIdentity();
      m.setRotationY(Math.cos(coeffs.phase(3, now) + x));
      m.multiply(makeRotX(Math.cos(coeffs.phase(2, now) + x)));
      m.multiply(makeTrans(0.0, 0.0, z * 0.45 * kScale));
      m.multiply(makeScale(sphereScale * 0.25, sphereScale * 0.25,
          sphereScale * (0.75 - 0.5 * z.abs())));
      s.ComputeInvTransform();

      n += _offset;
      x += dx;
    }
  }
}

class Orbit extends Form {
  final List<Torus> _tori = <Torus>[];

  Orbit() : super("Orbit") {
    for (int i = 0; i < 3; i++) {
      Torus t = new Torus(0.05 * kScale, (0.21 + i * 0.11) * kScale);
      addShape(t);
      _tori.add(t);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    {
      VM.Matrix4 m = _tori[0].transform;
      m.setIdentity();
      m.setRotationY(coeffs.value(8, now) * 10.0);
      m.multiply(makeRotX(coeffs.value(9, now) * 8.0));
      _tori[0].ComputeInvTransform();
    }

    {
      VM.Matrix4 m = _tori[1].transform;
      m.setIdentity();
      m.setRotationZ(coeffs.value(6, now) * 7.5);
      m.multiply(makeRotY(coeffs.value(7, now) * 6.0));
      _tori[1].ComputeInvTransform();
    }

    {
      VM.Matrix4 m = _tori[2].transform;
      m.setIdentity();
      m.setRotationX(coeffs.value(4, now) * 5.0);
      m.multiply(makeRotZ(coeffs.value(5, now) * 4.0));
      _tori[2].ComputeInvTransform();
    }
  }
}

class TorusBox extends Form {
  final List<Torus> _tori = <Torus>[];

  TorusBox() : super("TorusBox") {
    for (int i = 0; i < 3; i++) {
      Torus t = new Torus(0.05 * kScale, 0.4 * kScale);
      addShape(t);
      _tori.add(t);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    double a = coeffs.phase(1, now);
    {
      VM.Matrix4 m = _tori[0].transform;
      m.setIdentity();
      m.multiply(makeTrans(0.0, 0.0, 0.41 * kScale * Math.cos(a * 5.0)));
      _tori[0].ComputeInvTransform();
    }
    {
      VM.Matrix4 m = _tori[1].transform;
      m.setIdentity();
      m.setRotationX(HALFPI);
      m.multiply(makeTrans(
          0.0, 0.0, 0.41 * kScale * Math.cos((a + Math.PI * 0.3333) * 5.0)));
      _tori[1].ComputeInvTransform();
    }
    {
      VM.Matrix4 m = _tori[2].transform;
      m.setIdentity();
      m.setRotationY(HALFPI);
      m.multiply(makeTrans(
          0.0, 0.0, 0.41 * kScale * Math.cos((a + Math.PI * 0.6666) * 5.0)));
      _tori[2].ComputeInvTransform();
    }
  }
}

class RingOfTori extends Form {
  List<Torus> _tori = [];
  // >= 2
  RingOfTori(int n) : super("RingOfTori(${n})") {
    for (int i = 0; i < n; i++) {
      Torus t = new Torus(0.04 * kScale, 0.18 * kScale);
      _tori.add(t);
      addShape(t);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    double x = 0.0;
    final double dx = TWOPI / _shapes.length;
    for (Torus torus in _tori) {
      // _tori[i].SetScale3(0.75 + 0.25 * Math.cos(coeffs.phase(6, now) * 3.0), 1.0, 0.85);
      VM.Matrix4 m = torus.transform;
      m.setIdentity();
      m.setRotationZ(x);
      m.multiply(makeRotX(coeffs.phase(9, now) * 3.0));
      m.multiply(
          makeTrans(0.26 * kScale * coeffs.value(8, now, 3.0), 0.0, 0.0));
      m.multiply(makeRotY(coeffs.phase(7, now) * 3.0));
      m.multiply(makeScale(0.75 + 0.25 * coeffs.value(6, now, 3.0), 1.0, 0.85));
      torus.ComputeInvTransform();
      x += dx;
    }
  }
}

class SpheresAndCapsules extends Form {
  List<SphereUnit> _spheres = [];
  List<Capsule> _capsules = [];

  final sphereScale = 0.05 * kScale;

  SpheresAndCapsules(int n) : super("SpheresAndCapsules(${n})") {
    for (int i = 0; i < n; i++) {
      Capsule c = new Capsule(0.03 * kScale, 1.0 * kScale);
      _capsules.add(c);
      addShape(c);
    }

    for (int i = 0; i < n; i++) {
      SphereUnit s = new SphereUnit();
      _spheres.add(s);
      addShape(s);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    final double angle = TWOPI / _spheres.length;

    double a = 0.0;
    for (SphereUnit s in _spheres) {
      // TODO: good effect: swap transforms betweeb spheres and capsules
      VM.Matrix4 m = s.transform;
      m.setIdentity();
      m.setRotationY(coeffs.value(2, now, 4.0));
      m.multiply(makeRotX(coeffs.value(1, now, 4.0)));
      m.multiply(makeRotZ(a));
      m.multiply(makeTrans(0.4 * kScale * coeffs.value(0, now, 5.0), 0.0, 0.0));
      m.multiply(makeScale(sphereScale, sphereScale, sphereScale));
      s.ComputeInvTransform();
      a += angle;
    }

    a = 0.0;
    for (Capsule cap in _capsules) {
      cap.SetScale2(1.0, 0.07 + 0.07 * coeffs.value(6, now, 5.0));
      VM.Matrix4 m = cap.transform;
      m.setIdentity();
      m.setRotationY(coeffs.value(7, now, 4.0));
      m.multiply(makeRotX(coeffs.value(6, now, 4.0)));
      m.multiply(makeRotZ(a));
      m.multiply(
          makeTrans(0.35 * kScale * coeffs.value(5, now, 5.0), 0.0, 0.0));
      m.multiply(makeRotY(coeffs.value(4, now, 4.0)));
      m.multiply(makeRotX(coeffs.value(3, now, 4.0)));
      cap.ComputeInvTransform();
      a += angle;
    }
  }
}

class CubesAndCapsules extends Form {
  final List<RoundedCubeUnit> _cubes = <RoundedCubeUnit>[];
  final List<Capsule> _capsules = <Capsule>[];
  final double _offset;
  final double cubeScale = 0.03 * kScale;

  CubesAndCapsules(int n, this._offset)
      : super("CubesAndCapsules(${n}, ${_offset})") {
    for (int i = 0; i < n; i++) {
      Capsule c = new Capsule(0.03 * kScale, 1.0 * kScale);
      _capsules.add(c);
      addShape(c);
    }
    for (int i = 0; i < n; i++) {
      RoundedCubeUnit s = new RoundedCubeUnit();
      _cubes.add(s);
      addShape(s);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    final double hoffset = 0.35 * kScale;

    double x = 0.0;

    for (RoundedCubeUnit cube in _cubes) {
      VM.Matrix4 m = cube.transform;
      m.setIdentity();
      m.setRotationY(coeffs.value(7, now, 4.0));
      m.multiply(makeRotX(coeffs.value(6, now, 4.0)));

      m.multiply(makeTrans(
          hoffset * Math.cos(coeffs.phase(8, now) * 2.0 + x),
          hoffset * Math.cos(coeffs.phase(9, now) * 3.0 + x),
          hoffset * Math.cos(coeffs.phase(10, now) * 3.0 + x)));

      // @LOCALMOD prevent cube from getting zero size
      m.multiply(makeScale(
          cubeScale * (1.0 + 0.7 * coeffs.value(11, now, 5.0)),
          cubeScale * (1.0 + 0.7 * coeffs.value(12, now, 5.0)),
          cubeScale * (1.0 + 0.7 * coeffs.value(13, now, 5.0))));

      x += _offset;
      cube.ComputeInvTransform();
    }

    final double coffset = 0.32 * kScale;
    x = 0.0;
    for (Capsule cap in _capsules) {
      cap.SetScale2(1.0, 0.07 + 0.07 * coeffs.value(19, now, 5.0));
      VM.Matrix4 m = cap.transform;
      m.setIdentity();
      m.setRotationY(coeffs.value(15, now, 4.0));
      m.multiply(makeRotX(coeffs.value(14, now, 4.0)));
      m.setTranslationRaw(
          coffset * Math.cos(coeffs.phase(16, now) * 3.0 + x),
          coffset * Math.cos(coeffs.phase(17, now) * 2.0 + x),
          coffset * Math.cos(coeffs.phase(18, now) * 2.0 + x));
      x += _offset;
      cap.ComputeInvTransform();
    }
  }
}

class Octahedron extends Form {
  final sphereScale = 0.05 * kScale;
  final List<Capsule> _capsules = <Capsule>[];
  final List<SphereUnit> _spheres = <SphereUnit>[];

  Octahedron() : super("Octahedron") {
    for (int i = 0; i < 12; i++) {
      Capsule c = new Capsule(0.03 * kScale, 0.25 * kScale);
      addShape(c);
      _capsules.add(c);
    }

    for (int i = 0; i < 3; i++) {
      SphereUnit s = new SphereUnit();
      addShape(s);
      _spheres.add(s);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.3 * kScale;
    for (int i = 0; i < 3; i++) {
      VM.Matrix4 m = _spheres[i].transform;
      setTrans(
          m,
          offset * coeffs.value(i * 2 + 1, now, 5.0),
          offset * coeffs.value(i * 2 + 2, now, 5.0),
          offset * coeffs.value(i * 2 + 3, now, 5.0));
      m.multiply(makeScale(sphereScale, sphereScale, sphereScale));
      _spheres[i].ComputeInvTransform();
    }

    final double osize = 0.23 * kScale;
    final double angle = Math.PI * 0.25;

    for (int i = 0; i < 4; i++) {
      VM.Matrix4 m = _capsules[i].transform;
      m.setIdentity();
      m.setRotationY((i == 1 || i == 2) ? angle : -angle);
      m.setTranslationRaw(
          (i & 1 == 1) ? osize : -osize, 0.0, (i & 2 == 2) ? osize : -osize);
      _capsules[i].ComputeInvTransform();
    }

    for (int i = 4; i < 8; i++) {
      VM.Matrix4 m = _capsules[i].transform;
      m.setIdentity();
      m.setRotationX((i == 4 || i == 7) ? angle : -angle);
      m.setTranslationRaw(
          0.0, (i & 1 == 1) ? osize : -osize, (i & 2 == 2) ? osize : -osize);
      _capsules[i].ComputeInvTransform();
    }

    for (int i = 8; i < 12; i++) {
      VM.Matrix4 m = _capsules[i].transform;
      m.setIdentity();
      m.setRotationY(HALFPI);
      m.multiply(makeRotX((i == 8 || i == 11) ? angle : -angle));
      m.setTranslationRaw(
          (i & 1 == 1) ? osize : -osize, (i & 2 == 2) ? osize : -osize, 0.0);
      _capsules[i].ComputeInvTransform();
    }
  }
}

class UFO extends Form {
  List<SphereUnit> _spheres = [];
  Torus _torus;
  // 7, 10
  final double sphereScale = 0.04 * kScale;

  UFO(int num) : super("UFO(${num})") {
    for (int i = 0; i < num; i++) {
      SphereUnit s = new SphereUnit();
      _spheres.add(s);
      addShape(s);
    }
    _torus = new Torus(0.04 * kScale, 0.43 * kScale);
    addShape(_torus);
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.35 * kScale;
    final double p = TWOPI / _spheres.length;

    {
      final double s = 0.75 + 0.25 * coeffs.value(13, now, 3.0);
      VM.Matrix4 m = _torus.transform;
      m.setIdentity();
      m.setRotationX(coeffs.value(17, now) * 5.0);
      m.multiply(makeRotY(coeffs.value(16, now) * 4.0));
      m.multiply(makeRotZ(coeffs.value(15, now) * 3.0));
      m.multiply(makeScale(s, s, 1.0));
      _torus.ComputeInvTransform();
    }
    for (int i = 0; i < _spheres.length; i++) {
      final double phase = p * i;
      double angle = coeffs.phase(6, now) * 10.0 + phase;
      double ripple0 = (Math.cos(angle) + 1.0) * 0.5;
      ripple0 *= ripple0;
      ripple0 *= ripple0;
      ripple0 *= ripple0;
      double ripple1 = (Math.cos(angle + (Math.PI * 0.25)) + 1.0) * 0.5;
      ripple1 *= ripple1;
      ripple1 *= ripple1;
      ripple1 *= ripple1;
      //
      VM.Matrix4 m = _spheres[i].transform;
      // m.setIdentity();
      //m = makeRot(e100, coeffs.value(8, now)) * m;
      //m = makeRot(e010, coeffs.value(9, now)) * m;
      //
      final double x = offset * (0.75 + 0.25 * coeffs.value(12, now, 3.0));
      m.setIdentity();
      m.setRotationZ(phase);
      m.multiply(makeTrans(x, 0.0, 0.0));
      m.multiply(makeRotZ(coeffs.value(11, now)));
      m.multiply(makeRotX(coeffs.value(10, now)));
      double xz = sphereScale * (1.0 + ripple1 - (ripple0 * 0.5));
      double y = sphereScale * (1.0 + ripple0 - (ripple1 * 0.5));
      m.multiply(makeScale(xz, y, xz));
      _spheres[i].ComputeInvTransform();
    }
  }
}

abstract class KubeCommon extends Form {
  List<RoundedCubeUnit> _hexahedrons = [];
  SphereUnit _sphere;
  List<double> shift;
  final double sphereScale = 0.05 * kScale;

  KubeCommon(String name, int n, int s) : super(name) {
    for (int i = 0; i < n; i++) {
      RoundedCubeUnit h = new RoundedCubeUnit();
      _hexahedrons.add(h);
      addShape(h);
    }
    _sphere = new SphereUnit();
    addShape(_sphere);

    shift = new List<double>(s);
  }

  void UpdateShift(Coeffs coeffs, double now) {
    for (int i = 0; i < shift.length; i++) {
      double s = coeffs.value(i, now) * 5.0;
      s = s.clamp(-1.0, 1.0);
      if (s < 0.0) {
        s = (s + 1.0) * (s + 1.0) - 1.0;
      } else {
        s = 1.0 - (1.0 - s) * (1.0 - s);
      }
      shift[i] = s;
    }
  }

  void AnimateSphere(Coeffs coeffs, double now) {
    final double offsetSphere = 0.35 * kScale;
    VM.Matrix4 m = _sphere.transform;
    setTrans(
        m,
        offsetSphere * Math.sin(coeffs.phase(10, now) * 4.0),
        offsetSphere * Math.sin(coeffs.phase(11, now) * 4.0),
        offsetSphere * Math.sin(coeffs.phase(12, now) * 4.0));
    // vary up to 20%
    double scale = sphereScale * (1.0 + 0.2 * coeffs.value(3, now));
    // scale 1.0 means original size
    m.multiply(makeScale(scale, scale, scale));
    _sphere.ComputeInvTransform();
  }

  void RotateCube(VM.Matrix4 m, double a, double b, double c) {
    m.multiply(makeRotZ(c * HALFPI));
    m.multiply(makeRotY(b * HALFPI));
    m.multiply(makeRotX(a * HALFPI));
  }
}

class Kube extends KubeCommon {
  static final double cubeScale = 0.06 * kScale;

  Kube() : super("Kube", 9, 10);

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.18 * kScale;

    UpdateShift(coeffs, now);
    AnimateSphere(coeffs, now);

    for (int i = 0; i < 8; i++) {
      VM.Matrix4 m = _hexahedrons[i].transform;
      m.setIdentity();
      RotateCube(m, shift[0], shift[1], shift[2]);
      m.multiply(makeTrans(
          i < 4 ? -offset * shift[3] : offset * shift[3],
          i ~/ 2 % 2 == 0 ? -offset * shift[4] : offset * shift[4],
          i % 2 == 0 ? -offset * shift[5] : offset * shift[5]));
      RotateCube(m, shift[6], shift[7], shift[8]);
      m.multiply(makeScale(cubeScale, cubeScale, cubeScale));
      _hexahedrons[i].ComputeInvTransform();
    }
    {
      VM.Matrix4 m = _hexahedrons[8].transform;
      m.setIdentity();
      m.multiply(makeScale(cubeScale, cubeScale, cubeScale));
      _hexahedrons[8].ComputeInvTransform();
    }
  }
}

class Kube2 extends KubeCommon {
  static final double cubeScale = 0.06 * kScale;

  Kube2() : super("Kube2", 6, 15);

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.23 * kScale;

    UpdateShift(coeffs, now);
    AnimateSphere(coeffs, now);

    for (int i = 0; i < 6; i++) {
      double x = 0.0;
      double y = 0.0;
      double z = 0.0;
      switch (i) {
        case 0:
          x = -offset * shift[6];
          break;
        case 1:
          x = offset * shift[6];
          break;
        case 2:
          y = -offset * shift[7];
          break;
        case 3:
          y = offset * shift[7];
          break;
        case 4:
          z = -offset * shift[8];
          break;
        case 5:
          z = offset * shift[8];
          break;
      }
      final int index = i & 0xe;
      VM.Matrix4 m = _hexahedrons[i].transform;
      m.setIdentity();
      RotateCube(m, shift[0 + index], shift[1 + index], shift[2 + index]);
      m.multiply(makeTrans(x, y, z));
      RotateCube(m, 0.0, shift[9 + index], shift[10 + index]);
      m.multiply(makeScale(cubeScale, cubeScale, cubeScale));
      _hexahedrons[i].ComputeInvTransform();
    }
  }
}

class Kube3 extends KubeCommon {
  static final double cubeScale = 0.06 * kScale;

  Kube3() : super("Kube3", 6, 21);

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.23 * kScale;

    UpdateShift(coeffs, now);
    AnimateSphere(coeffs, now);

    void t(IsoShape s, int a, int b, int c, int d, int e, int f, double x,
        double y, double z) {
      VM.Matrix4 m = s.transform;
      m.setIdentity();
      RotateCube(m, shift[a], shift[b], shift[c]);
      m.setTranslationRaw(x, y, z);
      RotateCube(m, shift[d], shift[e], shift[f]);
      m.multiply(makeScale(cubeScale, cubeScale, cubeScale));
      s.ComputeInvTransform();
    }

    for (int i = 0; i < 6; i++) {
      if (i < 2) {
        double x = offset * shift[3];
        double y = offset * shift[4];
        double z = offset * shift[5];
        if (i == 0) {
          x = 0.0;
        } else {
          y = 0.0;
        }
        t(_hexahedrons[i], 0, 1, 2, 6, 7, 8, x, y, z);
      } else if (i < 4) {
        double x = offset * shift[12];
        if (i == 2) x *= -1.0;
        double y = x;
        double z = 0.0;
        t(_hexahedrons[i], 9, 10, 11, 14, 15, 16, x, y, z);
      } else {
        double x = -offset * shift[3];
        double y = -offset * shift[4];
        double z = -offset * shift[5];
        if (i == 4) {
          x = 0.0;
        } else {
          y = 0.0;
        }
        t(_hexahedrons[i], 17, 18, 19, 6, 7, 20, x, y, z);
      }
    }
  }
}

class Kube4 extends KubeCommon {
  static final double cubeScale = 0.06 * kScale;

  Kube4() : super("Kube4", 6, 12);

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.36 * kScale;
    final double offset2 = 0.095 * kScale;
    UpdateShift(coeffs, now);
    AnimateSphere(coeffs, now);
    for (int i = 0; i < 6; i++) {
      double n;
      double nn = Math.cos(coeffs.phase(6, now) * 4.0 + i * TWOPI / 6.0);
      if (nn >= 0.0)
        n = 1.0 - (1.0 - nn) * (1.0 - nn);
      else
        n = (1.0 + nn) * (1.0 + nn) - 1.0;
      n *= offset;

      List<double> o = [
        n,
        n,
        offset2,
        -offset2,
        offset2,
        -offset2,
        n,
        n,
        offset2,
        -offset2,
        offset2,
        -offset2
      ];

      VM.Matrix4 m = _hexahedrons[i].transform;
      m.setIdentity();
      if (i == 0 || i == 3) {
        m.setRotationY(shift[i * 2 + 1] * HALFPI);
      } else {
        m.setRotationZ(shift[i * 2 + 1] * HALFPI);
      }

      if (i == 1 || i == 4) {
        m.multiply(makeRotY(shift[i * 2] * HALFPI));
      } else {
        m.multiply(makeRotX(shift[i * 2] * HALFPI));
      }
      m.multiply(makeScale(cubeScale, cubeScale, cubeScale));
      m.setTranslationRaw(o[i + 0], o[i + 2], o[i + 4]);
      _hexahedrons[i].ComputeInvTransform();
    }
  }
}

class Grinder extends KubeCommon {
  final double _xScale;
  final double _yzScale;

  // 4, 5,6
  Grinder(int n)
      : _xScale = 0.02 * kScale,
        // choose size so gear teeth mesh together comfortably
        _yzScale = 0.4 * kScale / Math.pow(n, 1.1),
        super("Grinder($n)", 2 * n, 1 /*unused*/);

  @override
  void Animate(Coeffs coeffs, double now) {
    AnimateSphere(coeffs, now);

    final int halfCount = _hexahedrons.length ~/ 2;
    final double t = coeffs.value(3, now) * 6.0;

    for (int i = 0; i < halfCount; i++) {
      final double phase = TWOPI / halfCount * i + t;
      VM.Matrix4 m = _hexahedrons[i].transform;
      m.setIdentity();
      m.setRotationZ(phase);
      m.multiply(make9and12(-Math.sin(phase), 0.38 * kScale));
      m.multiply(makeScale(_xScale, _yzScale, _yzScale));
      _hexahedrons[i].ComputeInvTransform();
    }
    for (int i = halfCount; i < _hexahedrons.length; i++) {
      final double phase = TWOPI / halfCount * (i - halfCount + 0.5) + t;

      VM.Matrix4 m = _hexahedrons[i].transform;
      m.setIdentity();
      m.setRotationY(HALFPI);
      m.multiply(makeRotZ(phase));
      m.multiply(make9and12(Math.sin(phase), 0.38 * kScale));
      m.multiply(makeScale(_xScale, _yzScale, _yzScale));
      _hexahedrons[i].ComputeInvTransform();
    }
  }
}

class Hexahedron extends Form {
  List<SphereUnit> _spheres = [];
  List<Capsule> _capsules = [];
  final double sphereScale = 0.05 * kScale;

  Hexahedron() : super("Cube") {
    for (int i = 0; i < 3; i++) {
      SphereUnit s = new SphereUnit();
      _spheres.add(s);
      addShape(s);
    }
    for (int i = 0; i < 12; i++) {
      Capsule c = new Capsule(0.03 * kScale, 0.2 * kScale);
      _capsules.add(c);
      addShape(c);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.3 * kScale;

    for (int i = 0; i < _spheres.length; i++) {
      VM.Matrix4 m = _spheres[i].transform;
      setTrans(
          m,
          offset * coeffs.value(i * 2 + 1, now, 5.0),
          offset * coeffs.value(i * 2 + 2, now, 5.0),
          offset * coeffs.value(i * 2 + 3, now, 5.0));
      m.multiply(makeScale(sphereScale, sphereScale, sphereScale));
      _spheres[i].ComputeInvTransform();
    }

    // iterated over cube edges
    final double cs = 0.25 * kScale; // cubesize
    for (int i = 0; i < _capsules.length; i++) {
      VM.Matrix4 m = _capsules[i].transform;
      m.setIdentity();

      if (i < 4) {
        m.setTranslationRaw(i & 1 == 0 ? cs : -cs, i & 2 == 0 ? cs : -cs, 0.0);
      } else if (i < 8) {
        m.setRotationX(HALFPI);
        m.setTranslationRaw(i & 1 == 0 ? cs : -cs, 0.0, i & 2 == 0 ? cs : -cs);
      } else {
        m.setRotationY(HALFPI);
        m.setTranslationRaw(0.0, i & 1 == 0 ? cs : -cs, i & 2 == 0 ? cs : -cs);
      }
      _capsules[i].ComputeInvTransform();
    }
  }
}

class Cylinder extends Form {
  List<SphereUnit> _spheres = [];
  List<Capsule> _capsules = [];
  List<Torus> _tori = [];

  final double sphereScale = 0.05 * kScale;
  Cylinder() : super("Cylinder") {
    for (int i = 0; i < 3; i++) {
      SphereUnit s = new SphereUnit();
      _spheres.add(s);
      addShape(s);
    }

    for (int i = 0; i < 2; i++) {
      Torus t = new Torus(0.04 * kScale, 0.36 * kScale);
      _tori.add(t);
      addShape(t);
    }

    for (int i = 0; i < 3; i++) {
      // @LOCALMOD: original was 0.29
      Capsule c = new Capsule(0.04 * kScale, 0.295 * kScale);
      _capsules.add(c);
      addShape(c);
    }
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.3 * kScale;

    for (int i = 0; i < _spheres.length; i++) {
      VM.Matrix4 m = _spheres[i].transform;
      setTrans(
          m,
          offset * coeffs.value(i * 2 + 1, now, 5.0),
          offset * coeffs.value(i * 2 + 2, now, 5.0),
          offset * coeffs.value(i * 2 + 3, now, 5.0));
      m.multiply(makeScale(sphereScale, sphereScale, sphereScale));
      _spheres[i].ComputeInvTransform();
    }
    setTrans(_tori[0].transform, 0.0, 0.0, -0.38 * kScale);
    _tori[0].ComputeInvTransform();
    setTrans(_tori[1].transform, 0.0, 0.0, 0.38 * kScale);
    _tori[1].ComputeInvTransform();

    final double a = Math.PI * coeffs.value(9, now);
    List<double> angles = [
      Math.PI + a,
      Math.PI * 0.66666 - 2.0 * a,
      2.0 * a - Math.PI * 0.66666
    ];

    for (int i = 0; i < _capsules.length; i++) {
      VM.Matrix4 m = _capsules[i].transform;
      m.setIdentity();
      m.setRotationZ(angles[i]);
      m.multiply(makeTrans(0.36 * kScale, 0.0, 0.0));
      _capsules[i].ComputeInvTransform();
    }
  }
}

class Rings extends Form {
  final double cubeScale = 0.06 * kScale;

  List<RoundedCubeUnit> _hexahedrons = [];
  Torus _torus;
  Rings(int n) : super("Rings($n)") {
    for (int i = 0; i < n; i++) {
      RoundedCubeUnit h = new RoundedCubeUnit();
      _hexahedrons.add(h);
      addShape(h);
    }
    _torus = new Torus(0.04 * kScale, 0.3 * kScale);
    addShape(_torus);
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    //final double offset = 0.3 * kScale;
    final double da = TWOPI / _hexahedrons.length;
    double angle = 0.0;
    for (RoundedCubeUnit cube in _hexahedrons) {
      VM.Matrix4 m = cube.transform;
      m.setIdentity();
      m.setRotationZ(coeffs.value(3, now) * 4.0);
      // @LOCALMOD: original had another Z - Rotation here
      m.multiply(makeRotY(angle));
      _r.setIdentity();
      _r[1] = 0.5 * coeffs.value(0, now, 5.0);
      _r[2] = 0.5 * coeffs.value(1, now, 5.0);
      _r[4] = 0.5 * coeffs.value(2, now, 5.0);
      _r[12] = 0.32 * kScale;
      m.multiply(_r);
      m.multiply(makeScale(cubeScale, cubeScale, cubeScale));
      cube.ComputeInvTransform();
      angle += da;
    }

    {
      VM.Matrix4 m = _torus.transform;
      m.setIdentity();
      m.setRotationY(coeffs.value(11, now) * 5.0);
      m.multiply(makeRotX(coeffs.value(10, now) * 5.0));
      // only vary the radius
      _torus.SetScale2(1.0, 1.0 + 0.1 * coeffs.value(12, now));
      _torus.ComputeInvTransform();
    }
  }
}

/*
class Tennis extends Form {
 List<RoundedCube> _hexahedrons = [];
 Sphere _sphere;
 VM.Vector3 _pos = new VM.Vector3.zero();
 VM.Vector3 _dir = new VM.Vector3.zero();

 Tennis() {
   Math.Random rng = new Math.Random();
   _pos.x = 1.0 + rng.nextDouble();
   _pos.y = rng.nextDouble();
   _pos.z = rng.nextDouble();
   for (int i = 0; i < 3; i++) {
     RoundedCube h = new RoundedCube(0.03 * kScale, 0.0 * kScale, 0.1 * kScale, 0.1 * kScale);
     _hexahedrons.add(h);
     addShape(h);
   }
   _sphere = new Sphere(0.05 * kScale);
   addShape(_sphere);

 }

}
 */

class Tetrahedron extends Form {
  List<SphereUnit> _spheres = [];
  List<Capsule> _capsules = [];

  Tetrahedron() : super("Tetrahedron") {
    for (int i = 0; i < 3; i++) {
      SphereUnit s = new SphereUnit();
      _spheres.add(s);
      addShape(s);
    }

    for (int i = 0; i < 6; i++) {
      Capsule c = new Capsule(0.03 * kScale, 0.32 * kScale);
      _capsules.add(c);
      addShape(c);
    }
  }

  final double radius_edge = 0.27 * kScale;

  // tetrahedon has 6 edges
  final List<List<double>> _params = [
    [0.0, -0.61548, 0.0],
    [0.0, -0.61548, Math.PI * 0.666666],
    [0.0, -0.61548, Math.PI * -0.666666],
    [Math.PI * 0.5, Math.PI * 0.166666, Math.PI],
    [Math.PI * 0.5, Math.PI * 0.166666, Math.PI * 0.33333333],
    [Math.PI * 0.5, Math.PI * 0.166666, Math.PI * -0.33333333],
  ];

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.3 * kScale;
    final double sphereScale = 0.05 * kScale;
    for (int i = 0; i < _spheres.length; i++) {
      VM.Matrix4 m = _spheres[i].transform;
      setTrans(
          m,
          offset * coeffs.value(i * 2 + 1, now, 5.0),
          offset * coeffs.value(i * 2 + 2, now, 5.0),
          offset * coeffs.value(i * 2 + 3, now, 5.0));
      m.multiply(makeScale(sphereScale, sphereScale, sphereScale));
      _spheres[i].ComputeInvTransform();
    }
    // Note, the capsule are not changing over time - so we do not have
    // to run this over and over again

    for (int i = 0; i < _capsules.length; i++) {
      final double a = _params[i][0];
      final double b = _params[i][1];
      final double c = _params[i][2];

      VM.Matrix4 m = _capsules[i].transform;
      m.setIdentity();
      if (c != 0.0) m.multiply(makeRotZ(c));
      if (b != 0.0) m.multiply(makeRotY(b));
      m.multiply(makeTrans(radius_edge, 0.0, 0.0));
      if (a != 0.0) m.multiply(makeRotX(a));
      _capsules[i].ComputeInvTransform();
    }
  }
}

class HeartForm extends Form {
  List<SphereUnit> _spheres = [];
  Heart _heart;

  final double sphereScale = 0.05 * kScale;

  HeartForm() : super("Heart") {
    for (int i = 0; i < 0; i++) {
      SphereUnit s = new SphereUnit();
      _spheres.add(s);
      addShape(s);
    }
    _heart = new Heart(0.3 * kScale);
    addShape(_heart);
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.3 * kScale;

    for (int i = 0; i < _spheres.length; i++) {
      VM.Matrix4 m = _spheres[i].transform;
      setTrans(
          m,
          offset * coeffs.value(i * 2 + 1, now, 5.0),
          offset * coeffs.value(i * 2 + 2, now, 5.0),
          offset * coeffs.value(i * 2 + 3, now, 5.0));
      m.multiply(makeScale(sphereScale, sphereScale, sphereScale));
      _spheres[i].ComputeInvTransform();
    }
    {
      VM.Matrix4 m = _heart.transform;
      m.setIdentity();
      m.setRotationY(coeffs.value(20, now, 4.0));
      m.multiply(makeRotX(-HALFPI));
      _heart.ComputeInvTransform();
    }
  }
}

class KnotAndSpheres extends Form {
  List<SphereUnit> _spheres = [];
  Knot _knot;

  final double sphereScale = 0.06 * kScale;

  // (2, 3, 3);
  // KnotAndSpheres(2, 5, 3)
  // KnotAndSpheres(3, 2, 2)
  // KnotAndSpheres(3, 4, 2)
  // KnotAndSpheres(3, 5, 2)
  KnotAndSpheres(int coils, int twists, int spheres)
      : super("KnotAndSpheres($coils,$twists,$spheres)") {
    for (int i = 0; i < spheres; i++) {
      SphereUnit s = new SphereUnit();
      _spheres.add(s);
      addShape(s);
    }

    // _knot = new Knot(kScale * 0.2, kScale * 0.1, kScale * 0.05, 3, 2);
    _knot =
        new Knot(kScale * 0.28, kScale * 0.14, kScale * 0.04, coils, twists);
    addShape(_knot);
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    final double offset = 0.35 * kScale;
    for (int i = 0; i < _spheres.length; i++) {
      VM.Matrix4 m = _spheres[i].transform;
      setTrans(
          m,
          offset * coeffs.value(i + 7, now, 5.0),
          offset * coeffs.value(i + 8, now, 5.0),
          offset * coeffs.value(i + 9, now, 5.0));
      m.multiply(makeScale(sphereScale, sphereScale, sphereScale));
      _spheres[i].ComputeInvTransform();
    }

    double fade0 = coeffs.value(6, now) * 0.5 + 0.5;
    double fade1 = fade0 * (coeffs.value(7, now, 5.0) * 0.5 + 0.5);
    _knot.SetScale2(1.0 + 0.1 * fade1 / 0.28, 1.0 - fade1);
  }
}

class KnotAndTorus extends Form {
  List<Sphere> _spheres = [];
  Knot _knot;
  Torus _torus;

  KnotAndTorus(int coils, int twists) : super("KnotAndTorus($coils,$twists)") {
    _torus = new Torus(0.04 * kScale, 0.42 * kScale);
    addShape(_torus);
    _knot =
        new Knot(kScale * 0.28, kScale * 0.14, kScale * 0.04, coils, twists);
    addShape(_knot);
  }

  @override
  void Animate(Coeffs coeffs, double now) {
    VM.Matrix4 m = _torus.transform;
    m.setIdentity();
    m.setRotationX(coeffs.value(17, now) * 5.0);
    m.multiply(makeRotY(coeffs.value(16, now) * 4.0));
    m.multiply(makeRotZ(coeffs.value(15, now) * 3.0));
    m.multiply(makeScale(0.75 + 0.25 * coeffs.value(13, now, 3.0),
        0.75 + 0.25 * coeffs.value(14, now, 4.0), 1.0));
    _torus.ComputeInvTransform();

    _knot.ComputeInvTransform();

    double fade0 = coeffs.value(6, now) * 0.5 + 0.5;
    double fade1 = fade0 * (coeffs.value(7, now, 5.0) * 0.5 + 0.5);
    _knot.SetScale2(1.0 + 0.1 * fade1 / 0.28, 1.0 - fade1);
  }
}

// ============================================================
// ============================================================

typedef double Evaluator(VM.Vector3 v);

abstract class Scaler {
  double evaluator(VM.Vector3 v);
}

class EmergeScaler implements Scaler {
  // between 0,1 -1
  double _scale;
  Evaluator _eval;

  EmergeScaler(this._scale, this._eval);

  @override
  double evaluator(VM.Vector3 v) {
    double val = _eval(v);

    double t = ((_scale - 1.0) * 1.5 + v.x / 1000.0) * 10.0;
    t = t * t * t;
    if (t < 0.0) {
      return t + val;
    }
    return val;
  }
}

class ShrinkScaler implements Scaler {
  double _scale;
  Evaluator _eval;

  ShrinkScaler(this._scale, this._eval);

  @override
  double evaluator(VM.Vector3 v) {
    if (_scale < 0.0) {
      return -1.0;
    }
    VM.Vector3 w = v * (1 / (_scale * _scale));
    return _eval(w);
  }
}
