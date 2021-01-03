%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GNU Octave that combines slices
% into a 3D plot. Press any button
% to start the animation.
% youtube.com/OscarVeliz
% @author Oscar Veliz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fractal = figure('Name','3D Newton Fractal');
scatter3([],[],[]);
axis equal;
hold on;
xlabel('x');
ylabel('y');
zlabel('i');
domain = [-2 2; -2 2];
range = [2 2; -2 -2]; 
zero = [0 0; 0 0]; 
xi = imread('xi-plane.png');
surf(domain,zero,range,'CData',xi,'FaceColor','texturemap');
iy = imread('iy-plane.png');
surf(zero,domain,range,'CData',iy,'FaceColor','texturemap');
xy = imread('xy-plane.png');
surf(domain,range,zero,'CData',xy,'FaceColor','texturemap');
waitforbuttonpress;
for i=1:360
    camorbit(1,0,'camera')
    drawnow
end