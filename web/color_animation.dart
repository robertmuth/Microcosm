import 'dart:html';
import 'dart:core';
import 'dart:math' as Math;
import 'color.dart';

CanvasElement canvas = querySelector("#area");
PreElement pre = querySelector("#colors");

Math.Random rng = new Math.Random(1);

AnimatedColor ac = new AnimatedColor(rng);

void animate(num now) {
  ac.Update(now / 100);
  window.requestAnimationFrame(animate);
  DateTime begin = new DateTime.now();
  CanvasRenderingContext2D c = canvas.context2D;
  num w = canvas.width;
  //num h = canvas.height;

  CanvasGradient g = c.createLinearGradient(0, 0, w, 0);
  g.addColorStop(0, ac.Color(0));
  g.addColorStop(0.5, ac.Color(1));
  g.addColorStop(1, ac.Color(2));
  c
    ..fillStyle = g
    ..fillRect(0, 0, w, canvas.height);
  DateTime end = new DateTime.now();
  int msec = end.difference(begin).inMilliseconds;
  pre.innerHtml = "${ac.Color(0)}  ${ac.Color(1)}  ${ac.Color(2)}";
  print ("canvas update took: ${msec}ms");
/*
  CanvasGradient g1 = c.createLinearGradient(0, 0, w / 2, 0);
  g1.addColorStop(0, 'red');
  g1.addColorStop(1.0, 'blue');
  CanvasGradient g2 = c.createLinearGradient(w/2, 0, w, 0);
  g2.addColorStop(0.0, 'blue');
  g2.addColorStop(1, 'green');

  c
    ..fillStyle = g1
    ..fillRect(0, 0, w / 2, h)
    ..fillStyle = g2
    ..fillRect(w / 2, 0, w / 2, h);
     */
}

void main() {
  window.requestAnimationFrame(animate);
}
