library pc_shaders;

import 'package:chronosgl/chronosgl.dart';

const String aVertex2 = "aVertex2";

const String vTexCord = "vTexCord";

const String uTexCordMix = "uTexCordMix";
const String uLightPos0 = "uLightPos0";
const String uLightDiffuse0 = "uLightDiffuse0";
const String uLightAmbient0 = "uLightAmbient0";
const String uLightSpecular0 = "uLightSpecular0";

const String uLightPos1 = "uLightPos1";
const String uLightDiffuse1 = "uLightDiffuse1";
const String uLightAmbient1 = "uLightAmbient1";
const String uLightSpecular1 = "uLightSpecular1";

const String uFogColor = "uFogColor";
const String uFogEnd = "uFogEnd";
const String uFogScale = "uFogScale";

const String uCoeffRate = "uCoeffRate";
const String uCoeffPhase = "uCoeffPhase";
const String uShapeInvTransforms = "uShapeInvTransforms";
const String uShapeParams = "uShapeParams";
const String uAnimatedColors = "uAnimatedColors";

const int kMaxShapesPerForm = 20;

void Setup() {
  IntroduceNewShaderVar(aVertex2, new ShaderVarDesc(VarTypeVec3, ""));

  IntroduceNewShaderVar(vTexCord, new ShaderVarDesc(VarTypeVec3, ""));

  IntroduceNewShaderVar(uLightPos0, new ShaderVarDesc(VarTypeVec3, ""));
  IntroduceNewShaderVar(uLightPos1, new ShaderVarDesc(VarTypeVec3, ""));
  IntroduceNewShaderVar(uLightDiffuse0, new ShaderVarDesc(VarTypeVec3, ""));
  IntroduceNewShaderVar(uLightDiffuse1, new ShaderVarDesc(VarTypeVec3, ""));
  IntroduceNewShaderVar(uLightAmbient0, new ShaderVarDesc(VarTypeVec3, ""));
  IntroduceNewShaderVar(uLightAmbient1, new ShaderVarDesc(VarTypeVec3, ""));
  IntroduceNewShaderVar(uLightSpecular0, new ShaderVarDesc(VarTypeVec3, ""));
  IntroduceNewShaderVar(uLightSpecular1, new ShaderVarDesc(VarTypeVec3, ""));
  IntroduceNewShaderVar(uTexCordMix, new ShaderVarDesc(VarTypeVec3, ""));
  IntroduceNewShaderVar(
      uAnimatedColors, new ShaderVarDesc(VarTypeVec3, "", arraySize: 3));

  IntroduceNewShaderVar(uFogColor, new ShaderVarDesc(VarTypeVec3, ""));
  IntroduceNewShaderVar(uFogScale, new ShaderVarDesc(VarTypeFloat, ""));
  IntroduceNewShaderVar(uFogEnd, new ShaderVarDesc(VarTypeFloat, ""));
  IntroduceNewShaderVar(
      uCoeffPhase, new ShaderVarDesc(VarTypeFloat, "", arraySize: 25));
  IntroduceNewShaderVar(
      uCoeffRate, new ShaderVarDesc(VarTypeFloat, "", arraySize: 25));
  IntroduceNewShaderVar(uShapeInvTransforms,
      new ShaderVarDesc(VarTypeMat4, "", arraySize: kMaxShapesPerForm));
  IntroduceNewShaderVar(uShapeParams,
      new ShaderVarDesc(VarTypeVec4, "", arraySize: kMaxShapesPerForm));
}

/*
const float PI = 3.14159265359;
const float TWOPI = 6.28318530718;

mat4 MakeTranslation(vec3 t) {
    return mat4(vec4(1.0, 0.0, 0.0, 0.0),
                vec4(0.0, 1.0, 0.0, 0.0),
                vec4(0.0, 0.0, 1.0, 0.0),
                vec4(t, 0.0));
}

mat4 MakeRotationX(float radians) {
    float c = cos(radians);
    float s = sin(radians);
    return mat4(vec4(1.0, 0.0, 0.0, 0.0),
                vec4(0.0, c, s, 0.0),
                vec4(0.0, -s, c, 0.0),
                vec4(0.0, 0.0, 0.0, 0.0));
}

mat4 MakeRotationY(float radians) {
    float c = cos(radians);
    float s = sin(radians);
    return mat4(vec4(c, 0.0, -s, 0.0),
                vec4(0.0, 1.0, 0.0, 0.0),
                vec4(s, 0.0, c, 0.0),
                vec4(0.0, 0.0, 0.0, 0.0));
}

mat4 MakeRotationZ(float radians) {
    float c = cos(radians);
    float s = sin(radians);
    return mat4(vec4(c, s, 0.0, 0.0),
                vec4(-s, c, 0.0, 0.0),
                vec4(0.0, 0.0, 1.0, 0.0),
                vec4(0.0, 0.0, 0.0, 0.0));
}

float phase(int i, float now) {
    return mod(${uCoeffPhase}[i] + ${uCoeffRate}[i] * now, TWOPI);
}

float value(int i, float now) {
    return cos(phase(i, now));
}

float value(int i, float now, float freq) {
    return cos(phase(i, now) * freq);
}
 */

const String Helper = """
float SDF1Sphere(vec3 pos, float radius2_inv) { 
    return dot(pos, pos) *  radius2_inv;
}

float SDF1SphereUnit(vec3 pos) { 
    return dot(pos, pos);
}

float SDF1Torus(vec3 pos, float orbit, float thickness2_inv) { 
    float t = length(pos.xy) - orbit;
    return (t * t + pos.z * pos.z) * thickness2_inv;
}

float SDF1Capsule(vec3 pos, float radius2_inv, float halflength) { 
   float z = abs(pos.z) - halflength;
   if (z < 0.0) z = 0.0;
   vec3 p2 = vec3(pos.xy, z);
   return dot(p2, p2) *  radius2_inv;
}

float SDF1Cube(vec3 pos, float rx2_inv, float ry2_inv, float rz2_inv) { 
   float dx = (pos.x * pos.x) * rx2_inv;
   float dy = (pos.y * pos.y) * ry2_inv;
   float dz = (pos.z * pos.z) * rz2_inv;
   return max(dx, max(dy, dz));
}

float SDF1RoundedCubeUnit(vec3 pos) { 
   float dx = abs(pos.x) - 1.0;
   float dy = abs(pos.y) - 1.0;
   float dz = abs(pos.z) - 1.0;
   vec3 d = vec3(dx >= 0.0 ? dx : 0.0, dy >= 0.0 ? dy : 0.0, dz >= 0.0 ? dz : 0.0);
   return dot(d, d) * 4.0;
}

float SDF1RoundedCubeUnit(vec3 pos, float rx, float ry, float rz, float e2) { 
   float dx = abs(pos.x) - rx;
   float dy = abs(pos.y) - ry;
   float dz = abs(pos.z) - rz;
   vec3 d = vec3(dx >= 0.0 ? dx : 0.0, dy >= 0.0 ? dy : 0.0, dz >= 0.0 ? dz : 0.0);
   return dot(d, d) * e2;
}

float SDF1Knot(vec3 pos, int coils, float twists_over_coils, float lat_offset, 
               float thickness2, float r1, float r2) { 
    float tmp = length(pos.xy) - r1;
    float lat = atan(pos.y, pos.x) * twists_over_coils;
    float r = 0.0;
    for (int i = 0; i < coils; i++) {
      float lon = lat + lat_offset * float(i);
      float hor = tmp - cos(lon) * r2;
      float ver = pos.z - sin(lon) * r2;
      float d = hor * hor + ver * ver;
      if (d == 0.0) {
        return 0.0;
      }
      r += thickness2 / d;
    }
    if (r == 0.0) return 1000.0;
    return 1.0 / r;
}

float SDF1Heart(vec3 pos, float r) {
    float x = pos.x / r * 2.0;
    float y = pos.y / r;
    float z = pos.z / r;
    float xx = x * x;
    float yy = y * y;
    float zz = z * z;
    float zzz = zz * z;

    float a = 2.0 * xx + yy + zz - 1.0;
    a = a * a * a;
    float b = 0.1 * xx * zzz + yy * zzz;
    return (a - b) + 1.0;
}

const float epsilon = 0.05;

// forward declaration
float FormEval(vec3 pos);

vec3 FormNormal(vec3 pos) {
    float v0 = FormEval(pos);
    float vx = FormEval(vec3(pos.x - epsilon, pos.y, pos.z));
    float vy = FormEval(vec3(pos.x, pos.y - epsilon, pos.z));
    float vz = FormEval(vec3(pos.x, pos.y, pos.z - epsilon));
    return vec3(vx - v0, vy - v0, vz - v0);
}

/*
vec3 NoInterpolation(vec3 a, vec3 b) {
    return a + b;
}

vec3 Interpolate(vec3 p1, vec3 p2) {
    float v1 = FormEval(p1);
    float v2 = FormEval(p2);
    if (v1 > 0.0 || v2 < 0.0) {
       {
          vec3 t = p1;
          p1 = p2;
          p2 = t;
       }
       {
           float t = v1;
           v1 = v2;
           v2 = t;
       }
    }    
    float d = v1 / (v1 - v2);
    if (d < 0.0) d = 0.0;
    if (d > 1.0) d = 1.0;
    return mix(p1, p2, d);
}
*/
""";

// Disable Fog via
//glFogf(GL_FOG_START, 100.0f);
//  glFogf(GL_FOG_END, 1000.0f);

const String ShaderV = """

vec3 MakeColorMix(vec3 cordmix, vec3 pos, vec3 v0, vec3 n) {
     float x1 = length(v0)  * 0.001;
    float x2 = (v0.x + v0.y + v0.z) * 0.001 + 2.0;
    float y1 = dot(normalize(-pos), n) * 0.5 + 0.49;
    float y2 = n.x * 0.49 * 0.001 + 0.5;
    return vec3(clamp(mix(x1, x2, cordmix.x), 0.0, 1.0),
                clamp(mix(y1, y2, cordmix.y), 0.0, 1.0),
                cordmix.z);
}

void main(void) {    
    vec3 v0 = ${aPosition};
    // vec3 v0 =  NoInterpolation(${aPosition}, ${aVertex2});
    // vec3 v0 = Interpolate(${aPosition}, ${aVertex2});
    
    vec3 v = v0 * ${iaScale} + ${iaTranslation};
    
    /* normalizing likely not necessary */
    //vec3 n = normalize( (vec4(${aNormal}, 0.0)).xyz);
    
    vec3 n = normalize( (vec4(FormNormal(v0), 0.0)).xyz);
    ${vNormal} = n;
     
    /* ${vNormal} = ${uNormalMatrix} * ${aNormal}; */
    vec4 pos = vec4(v, 1.0);
    ${vPosition} = pos.xyz;
    gl_Position = ${uPerspectiveViewMatrix} * pos;
    /* vsEyeVec = -vsPosition.xyz; */
    ${vTexCord} = MakeColorMix( ${uTexCordMix}, pos.xyz, v0, n);
}
""";

const String ShaderF = """
vec4 diffuse;
vec3 specular;

void directionalLight(vec3 lightvec, 
                      vec3 cDiffuse, 
                      vec3 cAmbient, 
                      vec3 cSpecular, 
                      vec3 normal,
                      vec3 eyevec){
  // diffuse
  float norm_dot_dir = max(0.0, dot(normal, lightvec));
  diffuse.rgb += max(cDiffuse * norm_dot_dir, cAmbient);
  if(norm_dot_dir > 0.0){
    // phong specular
    vec3 reflectvec = reflect(-lightvec, normal);
    specular += cSpecular * pow(max(dot(eyevec, reflectvec), 0.0), 50.0);
  }
}

vec4 getColorFromTexture(vec3 v) {
    vec4 col1;
    vec4 col2;
    if (v.x <= 0.5) {
       col1.rgb = mix(${uAnimatedColors}[0], ${uAnimatedColors}[1], v.x * 2.0);    
    } else {
        col1.rgb = mix(${uAnimatedColors}[1], ${uAnimatedColors}[2], (v.x - 0.5) * 2.0);    
    }
    
    if (v.y <= 0.5) {
       col2.rgb = mix(${uAnimatedColors}[0], ${uAnimatedColors}[1], v.y * 2.0);    
    } else {
        col2.rgb = mix(${uAnimatedColors}[1], ${uAnimatedColors}[2], (v.y - 0.5) * 2.0);    
    }
    
    return mix(col1, col2, v.z);
}

void main(void) {
  diffuse = vec4(0.0, 0.0, 0.0, 1.0);
  specular = vec3(0.0, 0.0, 0.0);

  vec3 norm = normalize(vNormal);
  vec3 eye = normalize(-${vPosition});

  directionalLight(${uLightPos0}, 
                   ${uLightDiffuse0},
                   ${uLightAmbient0},
                   ${uLightSpecular0},
                   norm, eye);

  directionalLight(${uLightPos1}, 
                   ${uLightDiffuse1},
                   ${uLightAmbient1},
                   ${uLightSpecular1},
                   norm, eye);

  //diffuse = vec4(1.0, 1.0, 1.0, 1.0);
  diffuse *= getColorFromTexture(${vTexCord});
  //diffuse *= col1;
  // pre-multiplied alpha
  // diffuse.rgb *= diffuse.a;
  // smart specular addition
  diffuse.rgb += specular * (vec3(1.0, 1.0, 1.0) - diffuse.rgb);
  // Fog scale = 1 /(uFogEnd - <FogStart>)
  float f =  clamp((uFogEnd - length(${vPosition})) * uFogScale, 0.0, 1.0);
  diffuse = mix(vec4(uFogColor, 1.0), diffuse, f);
  ${oFragColor} = diffuse;
}

""";

final ShaderObject mcVertexShader = new ShaderObject("mcV")
  ..AddAttributeVars([aPosition, aVertex2, aNormal, iaTranslation, iaScale])
  ..AddVaryingVars([vNormal, vPosition, vTexCord])
  ..AddUniformVars([
    uPerspectiveViewMatrix,
    uTexCordMix,
    uShapeInvTransforms,
    uShapeParams,
  ])
  ..SetBody([Helper, ShaderV]);

final ShaderObject mcFragmentShader = new ShaderObject("mcF")
  ..AddVaryingVars([vNormal, vPosition, vTexCord])
  ..AddUniformVars([
    uLightPos0,
    uLightDiffuse0,
    uLightAmbient0,
    uLightSpecular0,
    uLightPos1,
    uLightDiffuse1,
    uLightAmbient1,
    uLightSpecular1,
    uAnimatedColors,
    uFogColor,
    uFogEnd,
    uFogScale
  ])
  ..SetBody([ShaderF]);
