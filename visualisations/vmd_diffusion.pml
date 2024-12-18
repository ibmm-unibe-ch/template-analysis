align 1pak, 2ubq;
align 1alf, 2ubq;

show cartoon, 1alf;
set_color blue, [0.9921568627450981, 0.8784313725490196, 0.42745098039215684];
color blue, 1alf;

show cartoon, 1pak;
set_color orange, [1, 0.6274509803921569, 0.2784313725490196];
color orange, 1pak;

show cartoon, 2ubq;
set_color yellow, [0.21568627450980393, 0.9058823529411765, 0.5137254901960784];
color yellow, 2ubq;

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
save different.png