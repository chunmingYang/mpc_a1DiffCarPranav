function animation(t,z,parms)

R = parms.R;
phi = linspace(0,2*pi);
x_circle = R*cos(phi);
y_circle = R*sin(phi);

if (parms.writeMovie)
    mov = VideoWriter(parms.nameMovie); 
    open(mov);
end

n = length(t);
for i = 1:n
    x_c = z(i,1); y_c = z(i,2); theta = z(i,3);
    [xe,ye] = ptC_to_ptP(x_c,y_c,theta,parms);
    ze(i,:) = [xe ye];
end

for i = 1:n
              
        theta=z(i,3); 
        x_robot = ze(i,1) + x_circle;
        y_robot = ze(i,2) + y_circle;
        x_dir = ze(i,1) + [0 R*cos(theta)];
        y_dir = ze(i,2) + [0 R*sin(theta)];

         
         plot(ze(1:i,1),ze(1:i,2),'r','Linewidth',2); hold on;
         light_blue = [176,224,230]/255; %from https://www.rapidtables.com/web/color/blue-color.html
         h1= patch(x_robot,y_robot,light_blue);  
         h2 = line(x_dir,y_dir,'Color','black','Linewidth',2);
           
          axis('equal');
          span = max([-min(min(ze(:,1:2))) max(max(ze(:,1:2)))]);
          span = max([span 2]);
          axis([-span span -span span]);
          grid on;

           pause(parms.delay)
           delete(h1);
           delete(h2);
           
           if (parms.writeMovie)
            axis off %does not show axis
            set(gcf,'Color',[1,1,1]) %set background to white
            writeVideo(mov,getframe);
           end
           
end % end for loop

if (parms.writeMovie)
    close(mov);
end


function [x_p,y_p] = ptC_to_ptP(x_c,y_c,theta,parms)
c = cos(theta); s = sin(theta); 
p = [parms.px; parms.py];
Xp = [c -s; s c]*p + [x_c; y_c];
x_p = Xp(1); y_p = Xp(2);


