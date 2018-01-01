import 'dart:math' as Math;
import 'package:vector_math/vector_math.dart';
import 'dart:typed_data';

class Coeffs {
  final int num;
  Float32List rate;
  Float32List phase;
  Float32List _phase_orig;
  Float32List value;

  Coeffs(Math.Random rng, this.num, double speed) {
    rate = new Float32List(num);
    _phase_orig = new Float32List(num);
    phase = new Float32List(num);
    phase = new Float32List(num);
    value = new Float32List(num);

    for (int i = 0; i < num; i++) {
      rate[i] = rng.nextDouble() * 0.015 * speed;
      _phase_orig[i] = rng.nextDouble() * 2 * Math.PI - Math.PI;
      phase[i] = 0.0;
      value[i] = 0.0;
    }
  }

  void Update(double now) {
    for (int i = 0; i < num; i++) {
      phase[i] = (_phase_orig[i] + rate[i] * now).remainder(2.0 * Math.PI);
      value[i] = Math.cos(phase[i]);
    }
  }
}

void main() {
  final Vector3 zero = new Vector3.zero();
  final Vector3 xyz = new Vector3(3.0, 5.0, 7.0);
  final Vector3 e1 = new Vector3(1.0, 0.0, 0.0);
  final Vector3 e2 = new Vector3(0.0, 1.0, 0.0);
  //final Vector3 e3 = new Vector3(0.0, 0.0, 1.0);
  final Matrix4 mId = new Matrix4.identity();
  {
    Vector3 v = mId * e1;
    if (v == e1) {
      assert(v == e1);
    }
  }

  {
    Matrix4 m = new Matrix4.identity();
    m.setTranslation(e1);
    assert(e1 == m * zero);
    assert((e1 + e2) == m * e2);
    m.setTranslation(e2);
    assert((e2 * 2.0) == m * e2);
  }

  {
    Matrix4 m = new Matrix4.identity();
    m.setTranslation(e1);
    Matrix4 i = new Matrix4.zero();
    i.copyInverse(m);
    print("m:\n$m");
    print("i:\n$i");
    print("product:\n${i*m}");
    Vector3 v = m * xyz;
    assert(i * v == xyz);
  }

  {
    const double kScale = 300.0;
    Math.Random rng = new Math.Random(1);
    Coeffs coeffs = new Coeffs(rng, 10, 10.0);
    coeffs.Update(1000.0);
    double x = 0.3;
    Quaternion q = new Quaternion.identity();
    Vector3 v = new Vector3.zero();
    Matrix4 r = new Matrix4.zero();
    Vector3 z = new Vector3.zero();
    Matrix4 m = new Matrix4.zero();
    m.setIdentity();
    v.setValues(
        0.0, 0.27 * kScale * coeffs.value[7], 0.27 * kScale * coeffs.value[8]);
    m.setTranslation(v);
    print("AFTER TRANS\n $m");
    q.setValues(coeffs.phase[9] * 2.0 + x, 1.0, 0.0, 0.0);
    r.setFromTranslationRotation(z, q);
    m = r * m;
    print("AFTER ROT\n $m");
  }
}
