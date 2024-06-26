hide everything, 2ubq

#
select first, /2ubq///1-100;
select end, /2ubq///100-153;

show cartoon, first_;
set_color blue, [0.2901960784313726, 0.3137254901960784, 0.592156862745098];
color blue, first_;

show cartoon, end;
set_color orange, [0.9529411764705882, 0.5411764705882353, 0.3686274509803922];
color orange, end;

set line_width, 3; # set line width
set ray_opaque_background, off; # turn on transparent alpha channel
set cartoon_dumbbell_width, 0.1; # default is 0.1 --> The helix horizontal length
set cartoon_dumbbell_length, 1.6; # default is 1.6 --> Control the vertical length of each helix
set cartoon_dumbbell_radius, 0.17; # default radius is 0.15
set cartoon_oval_quality, 1000; # default setting is 10
set cartoon_fancy_helices, 1; # Activate fancy helix representation
set ray_trace_mode, 1;  #3 -->  quantized color + black outline
set ray_trace_color, black; # set the colour of the outline line
set ray_trace_gain, 0.5; # set the width of the ouline line
set specular, off; # turn off all reflections/set opaque
set depth_cue, 0; # turn off depth cueing
bg_color white; # set background color to white 

ray 6512, 3333;
save coloured.png
