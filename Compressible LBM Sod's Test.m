clc;
clear ;
tic;
%initial conditions
rho(1:200,1:10) = 1;                                                    
rho(201:400,1:10)=0.125;
pre(1:200,1:10) = 1;
pre(201:400,1:10) = 0.1;
ux = zeros(400,10);
uy = zeros(400,10);

for t = 1:30 %number of iterations
f = zeros(400,10); %density distribution function
enef = zeros(400,10); %energy * density function
enea = zeros(400,10);
eneb = zeros(400,10);
xf = zeros(400,10);%x-component of velocity * density 
yf = zeros(400,10);%y component of velocity * density
for a = 1:400 % for each node on the lattice 
    for b = 1:10
        c = floor(ux(a,b)); %macroscopic transfer of particles
        d = floor(uy(a,b));
        x = [a+c a+c a+c+1 a+c+1];
        y = [b+d b+d+1 b+d+1 b+d];
        ratio = [(1-(ux(a,b)-c))*(1-(uy(a,b)-d)),(1-(ux(a,b)-c))*(uy(a,b)-d),(ux(a,b)-c)*(uy(a,b)-d),(ux(a,b)-c)*(1-(uy(a,b)-d))];
        ratio = ratio.*rho(a,b); %the array of fractions the density breaks into
        cs = sqrt(2*pre(a,b)/rho(a,b));
        c1 = floor(cs);
        c2=c1+1;
        d1 = 1/4 * ((c2*c2)-(cs*cs))/((c2*c2)-(c1*c1));
        d2 = 1/4 * ((cs*cs)-(c1*c1))/((c2*c2)-(c1*c1));
        pos = [c1 0 -c1 0 c2 0 -c2 0; 0 c1 0 -c1 0 c2 0 -c2];
        vdx = [ 0 0 1 1];
        vdy = [ 0 1 1 0 ];
        for p = 1:4 %microscopic streaming for all 4 nodes corresponding to the node in consideration
            xarr = [];
            yarr = [];
            reverse = [];
            for n=1:8 
                x1 = x(p) + pos(1,n);
                y1 = y(p) + pos(2,n);
                coor =  conv(x1,y1);
                xarr = [xarr coor(1)];
                yarr = [yarr coor(2)]; 
                reverse = [reverse coor(3)];
            end
            for i = 1:4
                f(xarr(i), yarr(i)) = f(xarr(i), yarr(i)) + (d1*ratio(p));
                %enef(xarr(i), yarr(i)) = enef(xarr(i), yarr(i)) + (0.5*d1*ratio(p)*((xarr(i)-a)^2+(xarr(i)-b)^2));
                enef(xarr(i), yarr(i)) = enef(xarr(i), yarr(i)) + (0.5*d1*ratio(p)*((pos(1,i)+c+vdx(p))^2+(pos(2,i)+d+vdy(p))^2))+(0.6*2.5*ratio(p)*d1*pre(a,b)/rho(a,b));
                enea(xarr(i), yarr(i)) = enea(xarr(i), yarr(i)) + (0.5*d1*ratio(p)*(c1^2));
                eneb(xarr(i),yarr(i)) = eneb(xarr(i),yarr(i)) + (0.5*d1*ratio(p)*(ux(a,b)^2));
                xf(xarr(i), yarr(i)) = xf(xarr(i), yarr(i)) + (d1*ratio(p)*(pos(1,i)+c+vdx(p))*reverse(i));
                yf(xarr(i), yarr(i)) = yf(xarr(i), yarr(i)) + (d1*ratio(p)*(pos(2,i)+d+vdy(p)));
            end
            for i = 5:8
               f(xarr(i), yarr(i)) = f(xarr(i), yarr(i)) + (d2*ratio(p));
               enea(xarr(i), yarr(i)) = enea(xarr(i), yarr(i)) + (0.5*d2*ratio(p)*(c2^2));
               eneb(xarr(i),yarr(i)) = eneb(xarr(i),yarr(i)) + (0.5*d1*ratio(p)*(ux(a,b)^2));
               %enef(xarr(i), yarr(i)) = enef(xarr(i), yarr(i)) + (0.5*ratio(p)*d2*((xarr(i)-a)^2+(xarr(i)-b)^2));
               enef(xarr(i), yarr(i)) = enef(xarr(i), yarr(i)) + (0.5*ratio(p)*d2*((pos(1,i)+c+vdx(p))^2+(pos(2,i)+d+vdy(p))^2))+(0.6*2.5*ratio(p)*d2*pre(a,b)/rho(a,b));
               xf(xarr(i), yarr(i)) = xf(xarr(i), yarr(i)) + (d2*ratio(p)*(pos(1,i)+c+vdx(p))*reverse(i));
               yf(xarr(i), yarr(i)) = yf(xarr(i), yarr(i)) + (d2*ratio(p)*(pos(2,i)+d+vdy(p)));
            end        
        end
    end
    
end
  %p=(gamma - 1)*rho*ene  gamma=2
rho = f;
ux = xf./f; 
uy = yf./f;
%ene2 = enea./f;
%ene = enet./f;
%internal = ene - (0.5*ux.*ux);
%ene2 = ene2/0.4;
%internal2 = ene2 - 0.5*ux.*(ux);

%internal2 = internal2/0.4;
ene = enef./f;
internal = ene - (0.5*ux);
pre = 0.4*internal.*f;


end
internal(1:100,1:10) = 2.5;
internal(300:400,1:10) = 2;
rho(1:100,1:10) = 1;
rho(300:400,1:10) = 0.125;
pre(1:100,1:10) = 1;
pre(300:400,1:10) = 0.1;
ux(1:100,1:10) = 0;
ux(300:400,1:10) = 0;

%internal3 = zeros(400,10);
%%for i = 1:200
  %  internal3(i,:) = internal2(i,:) + 1.5;
   
%end
%for i = 201:400
   % internal3(i,:) = internal2(i,:) +1.2;
%end
%plot(internal3(:,8));
















    
    

        