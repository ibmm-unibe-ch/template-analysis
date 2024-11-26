cealign 1org, 1rfd
set_color orange, [0.9607843137254902, 0.4745098039215686, 0.22745098039215686]; #https://venngage.com/blog/color-blind-friendly-palette/
set_color purple, [0.6627450980392157, 0.35294117647058826, 0.6313725490196078];
set_colo lightblue, [0.5215686274509804, 0.7529411764705882, 0.9764705882352941];
set_color darkblue, [0.058823529411764705, 0.12549019607843137, 0.5019607843137255];

color orange, 1org;
color darkblue, 1rfd;
color lightblue, 1fla;
color purple, 1noi;


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

set_view (\
     0.328686565,   -0.943073213,    0.050726570,\
     0.196892589,    0.120953828,    0.972937763,\
    -0.923685670,   -0.309801102,    0.225443915,\
     0.000000000,    0.000000000, -173.641296387,\
   223.475494385,  160.888427734,  261.796356201,\
   136.900115967,  210.382476807,  -20.000000000 )
   
ray 6512, 3333;
save pymol.png