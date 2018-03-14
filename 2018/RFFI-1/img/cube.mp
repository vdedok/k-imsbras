beginfig(0);


scalefactor := 10.000000 mm;
rot_psi := -60.000000;  % Rotation angle round z-axis (first rotation)
rot_phi := 15.000000;  % Rotation angle round y-axis (second rotation)
alpha := -24.146108;    % == arctan(sin(rot_phi)*tan(rot_psi))
beta  := -8.498781;    % == arctan(sin(rot_phi)/tan(rot_psi))

phi_source := 30.000000;
theta_source := 30.000000;
upper_value := 0.990000;
lower_value := 0.750000;
radius := scalefactor;
delta_rho := radius/50.000000;
delta_phi := 360.0/80.000000;


  path p;
  path equator;
  transform T;
  c1:=lower_value;
  c2:=upper_value-lower_value;
  nx_source := sind(theta_source)*cosd(phi_source);
  ny_source := sind(theta_source)*sind(phi_source);
  nz_source := cosd(theta_source);
  phistop := 360.0;
  rhostop := radius - delta_rho/2.0;

  for rho=0.0cm step delta_rho until rhostop:
    for phi=0.0 step delta_phi until phistop:
      rhomid := rho + delta_rho/2.0;
      phimid := phi + delta_phi/2.0;
      x1 := 70 + rho*cosd(phi);
      y1 := 70 + rho*sind(phi);
      x2 := 70 + (rho+delta_rho)*cosd(phi);
      y2 := 70 + (rho+delta_rho)*sind(phi);
      x3 := 70 + (rho+delta_rho)*cosd(phi+delta_phi);
      y3 := 70 + (rho+delta_rho)*sind(phi+delta_phi);
      x4 := 70 + rho*cosd(phi+delta_phi);
      y4 := 70 + rho*sind(phi+delta_phi);
      p:=makepath makepen ((x1,y1)--(x2,y2)--(x3,y3)--(x4,y4)--(x1,y1));
      quot := (rhomid/radius);
      nx_object := quot*cosd(phimid);
      ny_object := quot*sind(phimid);
      nz_object := sqrt(1-quot*quot);
      prod:=nx_object*nx_source+ny_object*ny_source
            +nz_object*nz_source;
      if prod < 0.0:
         value := c1;
      else:
         value := c1 + c2*prod*prod;
      fi
      fill p withcolor value[black,white];
    endfor
  endfor


lines = 5;
start_x = 10;
start_y = 10;
shift = 20;
dshift_xs = 10;
dshift_ys = 6;
dshift_x = lines * dshift_xs;
dshift_y = lines * dshift_ys;

side_len = 100;

end_x = start_x + side_len;
end_y = start_y + side_len;

for i = start_y step shift until end_y:
  draw (start_x, i)..(end_x, i);
endfor;

draw (start_x, start_y)..(end_x, start_y) withpen pencircle scaled 1.4;
draw (start_x, end_y)..(end_x, end_y) withpen pencircle scaled 1.4;

for i = start_x step shift until end_x:
  draw (i, start_y)..(i, end_y);
endfor;

draw (start_x, start_y)..(start_x, end_y) withpen pencircle scaled 1.4;
draw (end_x, start_y)..(end_x, end_y) withpen pencircle scaled 1.4;

for i = start_x step shift until end_x:
  draw (i, end_y)..(i + dshift_x, end_y + dshift_y);
endfor;
draw (start_x, end_y)..(start_x + dshift_x, end_y + dshift_y)  withpen pencircle scaled 1.4;
draw (end_x, end_y)..(end_x + dshift_x, end_y + dshift_y) withpen pencircle scaled 1.4;

for i = start_y step shift until end_y:
  draw (end_x, i)..(end_x + dshift_x, i + dshift_y);
endfor;
draw (end_x, start_y)..(end_x + dshift_x, start_y + dshift_y) withpen pencircle scaled 1.4;

x_min = end_x;
x_max = end_x + lines * dshift_xs;

for i = 1 step 1 until lines:
  draw (x_min + i * dshift_xs, start_y + i * dshift_ys)..(x_min + i * dshift_xs, start_y + i * dshift_ys + side_len);
endfor;
draw (x_min + lines * dshift_xs, start_y + lines * dshift_ys)..(x_min + lines * dshift_xs, start_y + lines * dshift_ys + side_len) withpen pencircle scaled 1.4;

for i = 1 step 1 until lines:
  draw (start_x + i * dshift_xs, start_y + i * dshift_ys + side_len)..(start_x + i * dshift_xs + side_len, start_y + i * dshift_ys + side_len);
endfor;
draw (start_x + lines * dshift_xs, start_y + lines * dshift_ys + side_len)..(start_x + lines * dshift_xs + side_len, start_y + lines * dshift_ys + side_len) withpen pencircle scaled 1.4;

draw (40, 40){dir 65}..{dir 0}(145, 100) dashed evenly withpen pencircle scaled 0.8;
pickup pencircle scaled 3pt;
draw (40,40);
draw (145,100);

%%label.bot("x", (40,40));
%%label.bot("y", (145,100));


endfig;

end