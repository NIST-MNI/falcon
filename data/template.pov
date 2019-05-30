/*#include "colors.inc"
#include "textures.inc"*/

/* ORIENTATIONS */
//rotate 360*clock*y
#declare ObjectRotationTop = <0,0,0>; 
#declare ObjectRotationLeft = <0,90,90>; 
#declare ObjectRotationBottom = <0,180,0>; 
#declare ObjectRotationRight = <0,-90,-90>; 
#declare ObjectRotationFront = <90,0,180>;
#declare ObjectRotationBack = <-90,0,0>;

#declare CameraDistance = 300;
#declare CameraLocation = <0,0,CameraDistance>;
#declare CameraRotation = <0,0,0>;
//#declare CameraLights = <0.7,0.7,0.7>;
#declare CameraLights = <0.9,0.9,0.9>;
#declare SceneLights = <0.5,0.5,0.5>;
//#declare SceneLights = <1,1,1>;

#declare LookAtTop = <0,0,0>;
#declare LookAtBottom = <0,0,0>;
#declare LookAtLeft = <0,0,0>;
#declare LookAtRight = <0,0,0>;
#declare LookAtFront = <0,0,0>;
#declare LookAtBack = <0,0,0>;

#declare SceneLightsOn=0;
#declare CameraLightOn=1;

global_settings { ambient_light rgb<0.9, 0.9, 0.9> }


camera {
  perspective

  right <-WIDTH/MIN,0,0>*0.7
  //right <-1,0,0>*0.7
  //right <-1.33,0,0>*0.7

  up <0,HEIGHT/MIN,0>*0.7
  //up <0,1,0>*0.7
  // Left right: -1.33 1
  // Top bottom: -1 1.33
  // Front back: -1 1.00

  location CameraLocation
  look_at LookAtANGLE
  //rotate CameraRotation
}

#if (CameraLightOn)

light_source {
  // Camera light
  CameraLocation
  color rgb CameraLights// shadowless
 // area_light <5, 0, 0>, <0, 0, 5>, 5, 5
 // adaptive 1
 // jitter

  rotate CameraRotation
}

#end

#if (SceneLightsOn)

light_source {
  // top light
  <0,200,0>
  color rgb SceneLights// shadowless
  rotate CameraRotation
}

light_source {
  // right light
  <200,0,0>
  color rgb SceneLights// shadowless
  //White
  rotate CameraRotation
}

light_source {
  // left light
  <-200,0,0>
  color rgb SceneLights// shadowless
  rotate CameraRotation
}

light_source {
  // bottom light
  <0,-200,0>
  color rgb SceneLights// shadowless
  rotate CameraRotation
}

light_source {
  // back light
  <0,0,200>
  color rgb SceneLights// shadowless
  rotate CameraRotation
}

light_source {
  // front light
  <0,0,-200>
  color rgb SceneLights// shadowless
  rotate CameraRotation
}

#end

//sky_sphere {
//    pigment {
//	color rgb <0,0,0.3>
//      }    
//  }
