#version 3.7;

sphere
{
	<0,0,0>, 5
	pigment
	{
		image_map
		{
//			jpeg "../img/sphere-2222-ABCBBACC-thumb.jpg"
			png "../sphere.png"
			map_type 1
			interpolate 2
		}
	}
	
	finish { phong 0.75 }
	
	rotate -y*clock*360
}

light_source
{
	<-66,66,-33> color rgb 1
}

camera
{
	right <1,0,0>
	up <0,1,0>
	location <0,5,-15>
	look_at <0,0,0>
	angle 40
}
