import 'dart:html' as HTML;
import 'dart:core';
import 'dart:math' as Math;
import 'dart:async';
import 'dart:web_gl' as WEBGL;
import 'dart:typed_data';

import 'shapes.dart' as Shape;

import 'package:vector_math/vector_math.dart' as VM;
import 'package:chronosgl/chronosgl.dart';
import 'option.dart';
import 'webutil.dart';
import 'color.dart';
import 'shaders.dart' as Shaders;
import 'logging.dart' as log;
import 'grid.dart' as GRID;

final bool _logTiminig = false;

double GetRandom(Math.Random rng, double a, double b) {
  return rng.nextDouble() * (b - a) + a;
}

final VM.Vector3 kColorWhite = new VM.Vector3(1.0, 1.0, 1.0);
final VM.Vector3 kColorBlack = new VM.Vector3(0.0, 0.0, 0.0);

double AngleDelta(double angleSrc, double angleDst) {
  double delta = angleDst - angleSrc;
  if (0.0 <= delta && delta <= 180.0) {
    return delta;
  }
  if (180.0 < delta) {
    return delta - 360.0;
  }
  // delta must be negative
  if (-180.0 <= delta) {
    return delta;
  }
  return delta + 360.0;
}

abstract class MyCamera {
  void animate(double elapsed, double now);
}

class RotatingCamera implements MyCamera {
  final Camera _camera;
  final double _cameraOrbit;
  final double _targetOrbit;
  final double _angleOffset;

  RotatingCamera(
      this._camera, this._cameraOrbit, this._targetOrbit, this._angleOffset);

  @override
  void animate(double elapsedMs, double nowMs) {
    double h = 0.0;
    double dir = nowMs;
    double cx = _cameraOrbit * Math.sin(dir);
    double cz = _cameraOrbit * Math.cos(dir);
    _camera.setPos(cx, h, cz);
    dir += _angleOffset;
    double tx = _targetOrbit * Math.sin(dir);
    double tz = _targetOrbit * Math.cos(dir);
    _camera.lookAt(new VM.Vector3(tx, 0.0, tz));
  }
}

class ManualCamera implements MyCamera {
  final OrbitCamera _orbit;
  final double _azimuthDelta;

  ManualCamera(this._orbit, this._azimuthDelta) {
    _orbit.polar = Math.PI / 8;
  }

  @override
  void animate(double elapsed, double now) {
    _orbit.azimuth += _azimuthDelta;
  }
}

final VM.Vector3 e100 = new VM.Vector3(1.0, 0.0, 0.0);
final VM.Vector3 e010 = new VM.Vector3(0.0, 1.0, 0.0);
final VM.Vector3 e001 = new VM.Vector3(0.0, 0.0, 1.0);

class RandomRotation {
  final List<double> _phase = [0.0, 0.0, 0.0];
  final List<double> _rate = [0.0, 0.0, 0.0];
  final List<double> _rot = [0.0, 0.0, 0.0];
  final VM.Matrix4 rotMat = new VM.Matrix4.identity();

  RandomRotation(Math.Random rng) {
    for (int i = 0; i < _phase.length; i++) {
      _phase[i] = rng.nextDouble() * Math.PI * 2.0;
      _rate[i] = 0.075 + 0.05 * rng.nextDouble();
    }
  }

  void animate(double elapsed, double now) {
    for (int i = 0; i < _phase.length; i++) {
      _phase[i] += _rate[i] * elapsed * 0.0001;
      _rot[i] += Math.sin(_phase[i]) * elapsed * 0.00012;
    }
  }

  void Rotate(VM.Matrix4 m) {
    m.rotate(e100, _rot[0]);
    m.rotate(e010, _rot[1]);
    m.rotate(e001, _rot[2]);
  }
}

class AutoCamera implements MyCamera {
  final Camera _camera;
  final RandomRotation _randRot;

  final VM.Vector3 _start = new VM.Vector3.zero();
  final VM.Vector3 _end = new VM.Vector3.zero();
  final Math.Random _rng;
  final VM.Matrix4 _rot = new VM.Matrix4.identity();

  //double _t = 1.0;

  AutoCamera(this._camera, Math.Random rng)
      : _rng = rng,
        _randRot = new RandomRotation(rng) {
    _randomPos(_end);
  }

  void _randomPos(VM.Vector3 v) {
    _end.x = _rng.nextDouble() * 400.0 - 200.0;
    _end.y = _rng.nextDouble() * 400.0 - 200.0;
    _end.z = -600.0 / 4.0;
  }

  @override
  void animate(double elapsed, double now) {
    _randRot.animate(elapsed, now);
    _rot.setIdentity();
    _randRot.Rotate(_rot);

    _start.x -= _rot[2] * 0.09 * elapsed;
    _start.y -= _rot[6] * 0.09 * elapsed;
    _start.z -= _rot[10] * 0.09 * elapsed;
    final double d = 1000.0;
    if (_start.x < -d) _start.x += 2.0 * d;
    if (_start.x > d) _start.x -= 2.0 * d;
    if (_start.y < -d) _start.y += 2.0 * d;
    if (_start.y > d) _start.y -= 2.0 * d;
    if (_start.z < -d) _start.z += 2.0 * d;
    if (_start.z > d) _start.z -= 2.0 * d;

    VM.Matrix4 m = _camera.transform;
    m.setIdentity();
    m[12] = _start.x;
    m[13] = _start.y;
    m[14] = _start.z;

    _randRot.Rotate(m);
  }
}

final HTML.Element gFps = HTML.querySelector("#fps");
Options gOptions;

void HandleCommand(String cmd, String param) {
  log.LogInfo("HandleCommand: ${cmd} ${param}");
  switch (cmd) {
    case "N":
      NextForm();
      break;
    case "A":
      Toggle(HTML.querySelector(".about"));
      break;
    case "C":
      Toggle(HTML.querySelector(".config"));
      gOptions.SaveToLocalStorage();
      //gPc.UpdateVisibility(gOptions);
      break;
    case "P":
      Toggle(HTML.querySelector(".performance"));
      break;
    case "R":
      gOptions.SaveToLocalStorage();
      HTML.window.location.hash = "";
      HTML.window.location.reload();
      break;
    case "A+":
      Show(HTML.querySelector(".about"));
      break;
    case "A-":
      Hide(HTML.querySelector(".about"));
      break;
    case "F":
      ToggleFullscreen();
      break;
    case "C-":
      Hide(HTML.querySelector(".config"));
      gOptions.SaveToLocalStorage();
      break;
    case "C+":
      Show(HTML.querySelector(".config"));
      break;
    case "X":
      String preset =
          (HTML.querySelector("#preset") as HTML.SelectElement).value;
      gOptions.SetNewSettings(preset);
      HTML.window.location.reload();
      break;
    default:
      break;
  }
}

List<Shape.Form> allForms = [
  new Shape.Metaballs(7),

  new Shape.KnotAndSpheres(3, 2, 3),
  new Shape.KnotAndSpheres(2, 3, 3),
  // new Shape.KnotAndSpheres(2,5,3),
  new Shape.KnotAndSpheres(3, 4, 2),
  //new Shape.KnotAndSpheres(3,5,2),
  new Shape.KnotAndSpheres(3, 5, 2),
  //
  new Shape.KnotAndTorus(2, 3),
  new Shape.KnotAndTorus(3, 2),
  new Shape.KnotAndTorus(3, 4),
  //
  new Shape.HeartForm(),
  new Shape.TriangleOfSpheres(5),
  new Shape.StringOfEllipsoids(8, 0.7),
  new Shape.Brain(6),
  new Shape.Flower(14, 1.0),
  new Shape.Orbit(),
  new Shape.TorusBox(),
// new Shape.KnotAndSpheres // not started
// new Shape.KnotAndTorus // not started
  new Shape.RingOfTori(5),
  new Shape.Tetrahedron(),
  new Shape.Hexahedron(),
  new Shape.Octahedron(),
  //new Shape.SpheresAndCapsules(5),
  new Shape.SpheresAndCapsules(6),
  new Shape.SpheresAndCapsules(7),
  new Shape.CubesAndCapsules(6, 0.6),
  new Shape.CubesAndCapsules(6, 1.0),
  new Shape.Kube(),
  new Shape.Kube2(),
  new Shape.Kube3(),
  new Shape.Kube4(),

  new Shape.Cylinder(),
  //new Shape.Rings(5),
  new Shape.Rings(6),
  new Shape.Rings(7),
  //new Shape.Grinder(4),
  new Shape.Grinder(5),
  new Shape.Grinder(6),
  new Shape.UFO(7),
  new Shape.UFO(10),
  //new Shape.Tennis(),
];

const String oRandom = "Random";
const String oDetail = "detail";
const String oFormDuration = "formDuration";
const String oColorSpeed = "colorSpeed";
const String oFormSpeed = "formSpeed";
const String oForm = "form";

const String oCameraMode = "cameraMode";
const String oHideAbout = "hideAbout";

const String oRandomSeed = "randomSeed";
const String oLogLevel = "logLevel";
const String oFov = "fov";
const String oMode = "mode";

void OptionsSetup() {
  // This must go first - so that when the options are initialized
  // The right form can be marked selected.
  HTML.SelectElement shapes = HTML.querySelector("#form");
  shapes.append(new HTML.OptionElement(data: oRandom, value: oRandom));
  for (Shape.Form f in allForms) {
    HTML.OptionElement o = new HTML.OptionElement(data: f.name, value: f.name);
    shapes.append(o);
  }

  gOptions = new Options("microcosm")
    ..AddOption(oDetail, "O", "20")
    ..AddOption(oFormDuration, "D", "30")
    ..AddOption(oColorSpeed, "D", "20.0")
    ..AddOption(oFormSpeed, "D", "20.0")
    ..AddOption(oForm, "O", "Random")
    ..AddOption(oCameraMode, "O", "manual")
    ..AddOption(oFov, "I", "50")
    ..AddOption(oMode, "O", "0")
    ..AddOption(oRandomSeed, "I", "0")
    ..AddOption(oHideAbout, "B", "false", true)
    // Only in debug mode
    ..AddOption(oLogLevel, "I", "0", true);

  gOptions.AddSetting("Standard", {
    oForm: "Random",
    oCameraMode: "manual",
    oFov: "50",
    oRandomSeed: "0",
    oMode: "0",
    oDetail: "15",
    oFormDuration: "40",
    oColorSpeed: "20",
    oFormSpeed: "20",
  });

  gOptions.AddSetting("Cylinder", {
    oForm: "Cylinder",
    oCameraMode: "rotate",
    oFov: "50",
    oRandomSeed: "0",
    oMode: "0",
    oDetail: "15",
    oFormDuration: "40",
    oColorSpeed: "20",
    oFormSpeed: "20",
  });

  gOptions.AddSetting("TorusBox", {
    oForm: "TorusBox",
    oCameraMode: "rotate",
    oFov: "50",
    oRandomSeed: "0",
    oMode: "0",
    oDetail: "15",
    oFormDuration: "40",
    oColorSpeed: "20",
    oFormSpeed: "20",
  });

  gOptions.AddSetting("KnotAndTorus", {
    oForm: "KnotAndTorus(2,3)",
    oCameraMode: "rotate",
    oFov: "50",
    oRandomSeed: "0",
    oMode: "0",
    oDetail: "15",
    oFormDuration: "40",
    oColorSpeed: "20",
    oFormSpeed: "20",
  });

  gOptions.AddSetting("SimpleKaleidoscope", {
    oForm: "Random",
    oCameraMode: "rotate",
    oFov: "50",
    oRandomSeed: "0",
    oMode: "1",
    oDetail: "15",
    oFormDuration: "40",
    oColorSpeed: "20",
    oFormSpeed: "20",
  });

  gOptions.AddSetting("ComplexKaleidoscope", {
    oForm: "Random",
    oCameraMode: "random",
    oFov: "50",
    oRandomSeed: "0",
    oMode: "3",
    oDetail: "12.5",
    oFormDuration: "40",
    oColorSpeed: "20",
    oFormSpeed: "20",
  });

  gOptions.ProcessUrlHash();

  log.gLogLevel = gOptions.GetInt(oLogLevel);
  log.gLogLevel = 0;

  HTML.SelectElement presets = HTML.querySelector("#preset");
  for (String name in gOptions.SettingsNames()) {
    HTML.OptionElement o = new HTML.OptionElement(data: name, value: name);
    presets.append(o);
  }

  if (gOptions.GetBool(oHideAbout)) {
    var delay = const Duration(seconds: 4);
    new Timer(delay, () => Hide(HTML.querySelector(".about")));
  }
}

double GetCurrentLevelOfDetail() {
  String s = gOptions.Get(oDetail);
  if (s == "") throw new Exception("bad lod ${s}");
  return double.parse(s);
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

void ExtractMeshDataUpdate(MeshData md, GRID.Grid grid, Shape.Form shape,
    GRID.DfsOption opt, double scale, bool computeNormals) {
  Shape.Scaler scaler = new Shape.EmergeScaler(scale, shape.evaluator);
  opt
    ..seqNo += 1
    ..eval = shape.evaluator;

  if (scale != 1.0) opt.eval = scaler.evaluator;

  List<int> seeds = [];
  for (VM.Vector3 v in shape.InteriorPoints()) {
    int s = grid.FindSurfacePointFromSeedBits(v, opt.eval);
    if (s != -1) seeds.add(s);
  }
  List<int> hot = [];
  grid.ValueUpdateDfs(seeds, true, hot, opt);

  GRID.MeshInfo info = grid.CubeMarch(hot, opt, computeNormals);

  DateTime begin = new DateTime.now();
  md.AddVertices(info.GetVertices());
  md.AddAttribute(Shaders.aVertex2, info.GetVertices2(), 3);

  md.AddFaces(info.GetFaces());
  if (computeNormals) {
    assert(info.nVertices3 == info.nNormals3);
    md.AddAttribute(aNormal, info.GetNormals(), 3);
  }
  DateTime end = new DateTime.now();
  int msec = end.difference(begin).inMilliseconds;
  if (_logTiminig) print("ActualMeshGeneration [${msec}ms]");
}

// ============================================================
//
// ============================================================
class FormRotator {
  final ChronosGL _cgl;
  final double _durationSec;
  final double _fadeInOutSec;
  final List<Shape.Form> _forms;
  final Map<String, RenderProgram> _programCache = <String, RenderProgram>{};

  RenderProgram currentProgram;
  Shape.Form currentForm;
  double currentScale;

  double _currentFormStartTime = 0.0;

  FormRotator(this._cgl, this._durationSec, this._fadeInOutSec, this._forms) {
    ChangeForm(_forms[0]);
    currentScale = 1.0;
  }

  void ChangeForm(Shape.Form form) {
    currentForm = form;
    if (!_programCache.containsKey(form.name)) {
      print("compile new shader for ${form.name}");
      final ShaderObject vertexShader = new ShaderObject("mcV")
        ..AddAttributeVars(
            [aPosition, Shaders.aVertex2, iaTranslation, iaScale])
        ..AddVaryingVars([vNormal, vPosition, Shaders.vTexCord])
        ..AddUniformVars([
          uPerspectiveViewMatrix,
          Shaders.uTexCordMix,
          Shaders.uShapeInvTransforms,
          Shaders.uShapeParams,
        ])
        ..SetBody(
            [Shaders.Helper, form.ExtractShaderCodeForEval(), Shaders.ShaderV]);

      _programCache[form.name] = new RenderProgram(
          form.name, _cgl, vertexShader, Shaders.mcFragmentShader);
    }
    currentProgram = _programCache[form.name];
  }

  void UpdateTime(double nowSec, Math.Random rng, String selected) {
    if (selected == oRandom) {
      final double pos = nowSec - _currentFormStartTime;
      if (_currentFormStartTime == 0.0 || pos >= _durationSec) {
        Shape.Form old = currentForm;
        Shape.Form f = currentForm;

        while (old == f) {
          f = _forms[rng.nextInt(_forms.length)];
        }
        ChangeForm(f);
        currentScale = 0.0;
        _currentFormStartTime = nowSec;
      } else if (pos < _fadeInOutSec) {
        currentScale = pos / _fadeInOutSec;
      } else if (pos >= _durationSec - _fadeInOutSec) {
        currentScale = (_durationSec - pos) / _fadeInOutSec;
      } else {
        currentScale = 1.0;
      }
    } else if (selected != currentForm.name) {
      currentForm = _forms[0];
      for (Shape.Form f in _forms) {
        if (f.name == selected) {
          ChangeForm(f);
          currentScale = 1.0;
          break;
        }
      }
    }
  }
}

void NextForm() {
  String formName = gOptions.Get(oForm);
  int i;
  for (i = 0; i < allForms.length; ++i) {
    if (allForms[i].name == formName) {
      break;
    }
  }
  i = (i + 1) % allForms.length;
  print("next form $i is: ${allForms[i].name}");
  gOptions.Set(oForm, allForms[i].name);
}

void AddInstancingData(MeshData md, int mode) {
  if (mode == 0) {
    final int count = 1;
    Float32List translations = new Float32List(count * 3);
    Float32List scalings = new Float32List(count * 3);

    translations[0] = 0.0;
    translations[1] = 0.0;
    translations[2] = 0.0;

    scalings[0] = 1.0;
    scalings[1] = 1.0;
    scalings[2] = 1.0;

    md.AddAttribute(iaScale, scalings, 3);
    md.AddAttribute(iaTranslation, translations, 3);
    return;
  }

  // intersting effect depending on overlap

  //final d = dimension / 2.0;
  final int k = mode - 1;
  final int count = (2 * k + 1) * (2 * k + 1) * (2 * k + 1) * 8;
  //print("@@@@count ${count}");
  Float32List translations = new Float32List(count * 3);
  Float32List scalings = new Float32List(count * 3);

  int pos = 0;
  final double spacing = 5000.0;
  final double offset = 500.0;
  for (int x = -k; x <= k; x++) {
    for (int y = -k; y <= k; y++) {
      for (int z = -k; z <= k; z++) {
        for (int i = 0; i < 8; i++) {
          double sx = i & 1 == 0 ? 1.0 : -1.0;
          double sy = i & 2 == 0 ? 1.0 : -1.0;
          double sz = i & 4 == 0 ? 1.0 : -1.0;
          scalings[pos + 0] = sx;
          scalings[pos + 1] = sy;
          scalings[pos + 2] = sz;
          translations[pos + 0] = x * spacing + offset * sx;
          translations[pos + 1] = y * spacing + offset * sy;
          translations[pos + 2] = z * spacing + offset * sz;
          pos += 3;
        }
      }
    }
  }
  assert(pos == count * 3);
  md.AddAttribute(iaScale, scalings, 3);
  md.AddAttribute(iaTranslation, translations, 3);
}

void main() {
  print("startup");
  if (!HasWebGLSupport()) {
    HTML.window.alert("Your browser does not support WebGL.");
    return;
  }
  OptionsSetup();
  Shaders.Setup();

  HTML.document.body.onKeyDown.listen((HTML.KeyboardEvent e) {
    log.LogInfo("key pressed ${e.keyCode} ${e.target.runtimeType}");
    if (e.target.runtimeType == HTML.InputElement) {
      return;
    }
    String cmd = new String.fromCharCodes([e.keyCode]);
    HandleCommand(cmd, "");
  });

  HTML.ElementList<HTML.Element> buttons =
      HTML.document.body.querySelectorAll("button");
  log.LogInfo("found ${buttons.length} buttons");

  buttons.onClick.listen((HTML.Event ev) {
    String cmd = (ev.target as HTML.Element).dataset['cmd'];
    String param = (ev.target as HTML.Element).dataset['param'];
    HandleCommand(cmd, param);
  });

  // querySelector("#area").onClick.listen((MouseEvent ev) {
  //   log.LogInfo("click area ${ev.target.runtimeType}");
  //   HandleCommand("C", "");
  // });

  final int rndSeed = gOptions.GetInt(oRandomSeed);
  final Math.Random rng = rndSeed == 0 ? new Math.Random() : new Math.Random(rndSeed);

  final HTML.CanvasElement canvas = HTML.querySelector("#area");
  final ChronosGL chronosGL = new ChronosGL(canvas);

  final AnimatedColor animatedColor = new AnimatedColor(rng);

  final FormRotator formRotator = new FormRotator(
      chronosGL, gOptions.GetDouble(oFormDuration), 3.0, allForms);

  final dimension = Shape.kScale / 2.0;

  print("Grid setup");
  final GRID.GridDimension gd = new GRID.GridDimension(
      new VM.Vector3(-dimension, -dimension, -dimension),
      new VM.Vector3(dimension, dimension, dimension),
      GetCurrentLevelOfDetail());
  GRID.Grid grid = new GRID.Grid(gd);
  //GRID.Grid grid2 = new GRID.Grid(lbf, utn, GetCurrentLevelOfDetail() * 2.0);
  GRID.DfsOption opt = new GRID.DfsOption()..seqNo = 666;
  Shape.Coeffs coeffs = new Shape.Coeffs(rng, 25, 10.0);

  final Material shapeMat = new Material("shapemat")
    ..SetUniform(Shaders.uLightPos0,
        new VM.Vector3(dimension, dimension, dimension)..normalize())
    ..SetUniform(Shaders.uLightPos1,
        new VM.Vector3(-dimension, -dimension, 0.0)..normalize())
    ..SetUniform(Shaders.uLightDiffuse0, kColorWhite)
    ..SetUniform(Shaders.uLightDiffuse1, new VM.Vector3(0.2, 0.2, 0.4))
    ..SetUniform(Shaders.uLightAmbient0, new VM.Vector3(0.2, 0.2, 0.2))
    ..SetUniform(Shaders.uLightAmbient1, kColorBlack)
    ..SetUniform(Shaders.uLightSpecular0, kColorWhite)
    ..SetUniform(Shaders.uLightSpecular1, new VM.Vector3(0.3, 0.3, 0.6))
    ..SetUniform(Shaders.uFogColor, kColorBlack)
    ..SetUniform(Shaders.uFogEnd, 20 * Shape.kScale)
    ..SetUniform(Shaders.uFogScale, 1 / (20 * Shape.kScale))
    // Not this gets updated under the hood via AnimatedColor::_UpdateMixer
    ..SetUniform(Shaders.uTexCordMix, animatedColor.GetMixer())
    ..SetUniform("uColor", new VM.Vector3(1.0, 1.0, 1.0));

  final MeshData md =
      formRotator.currentProgram.MakeMeshData("grid", WEBGL.TRIANGLES);

  ExtractMeshDataUpdate(md, grid, allForms[0], opt, 1.0, false);
  final int mode = int.parse(gOptions.Get(oMode));
  // we use mirroring in kaleidoscope mode where this does not work
  chronosGL.enable(GL_CULL_FACE);
  chronosGL.cullFace(GL_FRONT);

  AddInstancingData(md, mode);

  // Camera
  OrbitCamera orbit = new OrbitCamera(1500.0, 0.0, 0.0, canvas);
  orbit.mouseWheelFactor = -0.1;

  Perspective perspective = new Perspective(orbit, 10.0, 40000.0);

  Map<String, MyCamera> cameraAnimations = {
    "manual": new ManualCamera(orbit, 0.0),
    "rotate": new ManualCamera(orbit, 0.01),
    "random": new AutoCamera(orbit, rng),
    //
  };

  String lastCameraMode = "";

  void animateCamera(double elapsed, double now) {
    //String mode = gOptions.Get(oCameraMode);
    // FIXME
    String mode = gOptions.Get(oCameraMode);
    if (cameraAnimations.containsKey(mode)) {
      if (mode != lastCameraMode) {
        lastCameraMode = mode;
      }
      cameraAnimations[mode].animate(elapsed, now);
    } else {
      log.LogError("unknown camera mode ${mode}");
    }
    perspective.UpdateFov(0.0 + gOptions.GetInt(oFov));
  }

  Float32List shapeInvTransforms =
      new Float32List(16 * Shaders.kMaxShapesPerForm);
  Float32List shapeParams = new Float32List(4 * Shaders.kMaxShapesPerForm);

  log.LogInfo("Starting ChronosGL main loop");

  double _lastTimeMs = 0.0;

  void animate(num timeMs) {
    timeMs = 0.0 + timeMs;
    double elapsed = timeMs - _lastTimeMs;
    _lastTimeMs = timeMs;

    // Camera Stuff
    orbit.azimuth += 0.001;
    orbit.animate(elapsed);

    animateCamera(elapsed, timeMs);

    // Update colors
    final double colorSpeed = gOptions.GetDouble(oColorSpeed);
    // This will update animatedColor.GetMixer()
    animatedColor.Update(colorSpeed * timeMs * 0.001);
    shapeMat.ForceUniform(Shaders.uAnimatedColors, animatedColor.colors);

    // Update form shape
    formRotator.UpdateTime(timeMs * 0.001, rng, gOptions.Get(oForm));
    final Shape.Form shape = formRotator.currentForm;

    // Animated Shape and set some uniforms updated by this
    final double formSpeed = gOptions.GetDouble(oFormSpeed);
    shape.Animate(coeffs, timeMs * formSpeed * 0.0001);
    shape.UpdateShapeInvTransform(shapeInvTransforms);
    shape.UpdateShapeParams(shapeParams);
    shapeMat.ForceUniform(Shaders.uShapeInvTransforms, shapeInvTransforms);
    shapeMat.ForceUniform(Shaders.uShapeParams, shapeParams);

    // Update Meshdata
    DateTime begin = new DateTime.now();
    ExtractMeshDataUpdate(
        md, grid, shape, opt, formRotator.currentScale, false);
    DateTime end = new DateTime.now();
    int msec = end.difference(begin).inMilliseconds;
    if (_logTiminig) print("UpdateMeshData ${md} [${msec}ms]");

    List<DrawStats> stats = [];
    formRotator.currentProgram.Draw(md, [perspective, shapeMat], stats);
    List<String> out = [];
    for (DrawStats d in stats) {
      out.add(d.toString());
    }

    UpdateFrameCount(timeMs, gFps, out.join("\n"));

    HTML.window.animationFrame.then(animate);
  }

  void resolutionChange(HTML.Event ev) {
    canvas.width = HTML.window.innerWidth;
    canvas.height = HTML.window.innerHeight;
    int w = canvas.clientWidth;
    int h = canvas.clientHeight;
    canvas.width = w;
    canvas.height = h;
    print("size change $w $h");
    perspective.AdjustAspect(w, h);
    chronosGL.viewport(0, 0, w, h);
  }

  resolutionChange(null);
  HTML.window.onResize.listen(resolutionChange);

  animate(0.0);
}
