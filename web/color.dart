library color;

import 'dart:math' as Math;
import 'dart:typed_data';
import 'package:vector_math/vector_math.dart' as VM;

VM.Vector3 HslToRgb(double h, double s, double l) {
  double r, g, b;

  if (h < 1.0 / 6.0) {
    // full red, some green
    r = 1.0;
    g = h * 6.0;
    b = 0.0;
  } else {
    if (h < 0.5) {
      // full green
      g = 1.0;
      if (h < 1.0 / 3.0) {
        // some red
        r = 1.0 - ((h - 1.0 / 6.0) * 6.0);
        b = 0.0;
      } else {
        // some blue
        b = (h - 1.0 / 3.0) * 6.0;
        r = 0.0;
      }
    } else {
      if (h < 5.0 / 6.0) {
        // full blue
        b = 1.0;
        if (h < 2.0 / 3.0) {
          // some green
          g = 1.0 - ((h - 0.5) * 6.0);
          r = 0.0;
        } else {
          // some red
          r = (h - 0.666667) * 6.0;
          g = 0.0;
        }
      } else {
        // full red, some blue
        r = 1.0;
        b = 1.0 - ((h - 5.0 / 6.0) * 6.0);
        g = 0.0;
      }
    }
  }
  // saturation influence
  r = 1.0 - (s * (1.0 - r));
  g = 1.0 - (s * (1.0 - g));
  b = 1.0 - (s * (1.0 - b));

  // luminosity influence
  r *= l;
  g *= l;
  b *= l;
  return new VM.Vector3(r, g, b);
}

VM.Vector3 _MakeSkew(Math.Random rng) {
  return new VM.Vector3(rng.nextDouble() * Math.PI * 2.0,
      rng.nextDouble() * Math.PI * 2.0, rng.nextDouble() * Math.PI * 2.0);
}

VM.Vector3 AnimateTexture(double now, double speed, VM.Vector3 skew) {
  // oscillator for object space texcoord
  double a = skew.x + now * 0.0181 * speed;
  // oscillator for eye space texcoord
  double b = skew.y + now * 0.0221 * speed;
  // oscillator for moving between object and eye space texcoords
  double c = skew.z + now * 0.0037 * speed;

  double x = Math.sin(a) * 0.5 + 0.5;
  double y = Math.sin(b) * 0.6 + 0.5;
  double z = Math.sin(c) + 0.5;
  return new VM.Vector3(x.clamp(0, 1), y.clamp(0, 1), z.clamp(0, 1));
}

// ===============================================================

const int _kColors = 3;

class AnimatedColor {
  final VM.Vector3 _skew;
  final VM.Vector3 _mixer = new VM.Vector3.zero();
  final Uint8List _color = new Uint8List(_kColors * 3);
  final Float32List _coeff_value = new Float32List(_kColors * 6);
  final Float32List _coeff_rate = new Float32List(_kColors * 6);
  final Float32List _coeff_phase = new Float32List(_kColors * 6);

  final Float32List _colors = new Float32List(_kColors * 3);

  AnimatedColor(Math.Random rng) : _skew = _MakeSkew(rng) {
    for (int i = 0; i < _kColors * 6; i++) {
      _coeff_phase[i] = rng.nextDouble() * Math.PI * 2.0;
      _coeff_rate[i] = rng.nextDouble() * 0.002 + 0.002;
    }
  }

  String Color(int i) {
    return "rgb(${_color[3* i +0]}, ${_color[3* i +1]}, ${_color[3* i +2]})";
  }

  Float32List get colors => _colors;

  VM.Vector3 GetMixer() {
    return _mixer;
  }

  void _UpdateMixer(double now) {
    // oscillator for object space texcoord
    double a = _skew.x + now * 0.0181;
    // oscillator for eye space texcoord
    double b = _skew.y + now * 0.0221;
    // oscillator for moving between object and eye space texcoords
    double c = _skew.z + now * 0.0037;

    _mixer.x = Math.sin(a) * 0.5 + 0.5;
    _mixer.y = Math.sin(b) * 0.6 + 0.5;
    _mixer.z = Math.sin(c) + 0.5;
  }

  void _UpdateCoeffs(double now) {
    for (int i = 0; i < _kColors * 6; i++) {
      _coeff_value[i] =
          (_coeff_phase[i] + now * _coeff_rate[i]) % (2.0 * Math.PI);
    }
  }

  static double GetHue(int i, double m) {
    switch (i % 3) {
      case 0:
        return m * 2.0;
      case 1:
        return 2.0 - m * 2.0;
      case 2:
        return m * 3.0;
      default:
        assert(false);
        return 0.0;
    }
  }

  void _UpdateColor(int i) {
    int c = i * 6;
    final double hue = GetHue(i, _coeff_value[c + 0] / (Math.PI * 2.0));
    // Saturation and luminosity are bonehead simple.
    final double sat = _coeff_value[c + 1] * 0.5 + 0.5;
    //double lum = _coeff_value[c + 2] * 0.5 + 0.5;
    final double lum = _coeff_value[c + 2] * 0.3 + 0.7;

    // The resulting colors here are a lot of earth tones:
    // off whites, and pastels.  Ick

    // Convert to rgb color space
    VM.Vector3 rgb = HslToRgb(hue, sat, lum);

    // Randomize a little more (now the color is *really* random)
    // Set this constant too high and there will be too many pure primaries and secondaries.
    double rgb_variance = 0.5;
    rgb.r += Math.cos(_coeff_value[c + 3]) * rgb_variance;
    rgb.g += Math.cos(_coeff_value[c + 4]) * rgb_variance;
    rgb.b += Math.cos(_coeff_value[c + 5]) * rgb_variance;
    rgb.clampScalar(0.0, 1.0);

    _colors[3 * i + 0] = rgb.x;
    _colors[3 * i + 1] = rgb.y;
    _colors[3 * i + 2] = rgb.z;

    _color[3 * i + 0] = (rgb.x * 256.0).floor().clamp(0, 255);
    _color[3 * i + 1] = (rgb.y * 256.0).floor().clamp(0, 255);
    _color[3 * i + 2] = (rgb.z * 256.0).floor().clamp(0, 255);
  }

  void Update(double now) {
    _UpdateMixer(now);
    _UpdateCoeffs(now);
    for (int i = 0; i < _kColors; i++) {
      _UpdateColor(i);
    }
  }
}
