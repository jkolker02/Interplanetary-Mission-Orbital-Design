clear, clc, close all

%----------PART 1-----------
%
%
%
%
%
%

figure('Name', 'Part 1: Transfer Orbit to Mars', 'Position', [100 300 600 600], 'NumberTitle','off')

[posE, velE] = planetEphemeris(juliandate(2004, 6, 5, 1, 52, 21),'Sun','Earth', '405', 'km');
[posM, velM] = planetEphemeris(juliandate(2005, 5, 14, 13, 23, 33), 'Sun', 'Mars', '405', 'km');
[posMPLOTTING, velMPLOTTING] = planetEphemeris(juliandate(2004, 6, 5, 1, 52, 21),'Sun','Mars', '405', 'km');

rE = 6378; %Radius earth
rM = 3390; %Radius Mars
v1 = sqrt(398600/(4*rE)); %Velocity in earth parking orbit
v2 = sqrt(42828/(4*rM)); %Velocity in mars capture orbit


ToF = 343 + 11.5/24; %time between the launch and arrival dates in days


[v1D, v2A] = lambert(posE, posM, ToF, 0, 132712440018)


v1inf = v1D - velE; %v infinity departing earth
v2inf = v2A - velM; % v infinity arriving at mars
magv1inf = norm(v1inf); % magnitude of v1 infinity
magv2inf = norm(v2inf); % magnitude of v2 infinity

vp = sqrt(magv1inf^2 + 2*398600/(4*rE)); %Velcoity of spacecraft at perigee on hyperbola at launch
vp2 = sqrt(magv2inf^2 + 2*42828/(4*rM)); %Velocity of spacecraft at pergiee on hyperbola at arrival at mars


deltav1 = vp - v1; %dvel = velocity on hyperbola at perigee - velocity on circular orbit
deltav2 = vp2 - v2;
magdV1 = norm(deltav1);
magdV2 = norm(deltav2);



disp("--------Part 1--------")
disp("Magnitude of delta V1: " + magdV1 + "km/s")
disp("Magnitude of delta V2: " + magdV2 + "km/s")
disp("Velcity upon entering Mars SOI: " + magv2inf + "km/s")



%--------PLOTTING------

X = [posE,v1D]; %All the initial parameters in vector 'X'
time = 29676600; %Desired Time elapsed
time2 = 45000000;
timespan = [0,time]; %Start and end times for ode45 to use
ode2=@drdr2;
options = odeset('AbsTol', 1e-9, 'RelTol', 1e-9, 'MaxStep',80000); %Found to help keep tolerances from drifting values


[t,y]= ode45(ode2,timespan,X, options); %Runs the ODE calculation numerically using a reference to the drdr function, time span, and initial conditions

X2 = [posE, velE];
[t2,y2] = ode45(ode2, timespan, X2, options);

xr = y(:,1); %Grabs x values of r
yr = y(:,2); %Grabs y values of r
zr = y(:,3); %Grabs z values of r

plot3(xr,yr,zr, 'linewidth', 2.5, 'Color', 'Black'); %Plots all values of r on 3d graph
hold on;

xr2 = y2(:,1); %Grabs x values of r
yr2 = y2(:,2); %Grabs y values of r
zr2 = y2(:,3); %Grabs z values of r
plot3(xr2, yr2, zr2, 'linewidth', 2.5,'Color', 'Blue')

plot3(0,0,0, '*', 'MarkerSize', 18, 'LineWidth', 2 , 'Color', '#ffa500'); %Plots the sun

X3 = [posMPLOTTING, velMPLOTTING];
[t3,y3] = ode45(ode2, timespan, X3, options);

xr3 = y3(:,1); %Grabs x values of r
yr3 = y3(:,2); %Grabs y values of r
zr3 = y3(:,3); %Grabs z values of r
plot3(xr3, yr3, zr3, 'linewidth',2.5, 'Color', 'Red')

view(0,89.9)
axis square
xlim([-3 * 10^8, 3 * 10^8])
ylim([-3 * 10^8, 3 * 10^8])






%------------PART 2------------
%
%
%
%

figure('NumberTitle','off', 'Name', 'Part 2: Transfer orbit to Venus, Mars', 'Position', [750 300 600 600]);

rV = 6051.8;
%V1 remains the same, same circular parking orbit from part 1 around earth
%
%rP = 4 radius of venus circular orbit. 

v2_2 = sqrt(324859/(4*rV));

[posV, velV] = planetEphemeris(juliandate(2004, 11, 20, 15, 10, 59),'Sun','Venus', '405', 'km');
[posE, velE] = planetEphemeris(juliandate(2004, 6, 5, 1, 52, 21),'Sun','Earth', '405', 'km');
[posVPLOTTING, velVPLOTTING] = planetEphemeris(juliandate(2004, 6, 5, 1, 52, 21),'Sun','Venus', '405', 'km');


ToF2 = 167.445; %days to go from earth to venus
[v1D2, v2A2] = lambert(posE, posV, -ToF2, 0, 132712440018);



v1inf2 = v1D2 - velE; %v infinity departing earth
v2inf2 = v2A2 - velV; % v infinity arriving at venus
magv1inf2 = norm(v1inf2); % magnitude of v1 infinity
magv2inf2 = norm(v2inf2); % magnitude of v2 infinity

velp = sqrt(magv1inf2^2 + 2*398600/(4*rE)); %Velcoity of spacecraft at perigee on hyperbola at launch from earth
velp2 = sqrt(magv2inf2^2 + 2*324859/(4*rV)); %Velocity of spacecraft at pergiee on hyperbola at arrival at mars


deltav1_2 = velp - v1; %dvel = velocity on hyperbola at perigee - velocity on circular orbit.  Delta V at earth
deltav2_2 = velp2 - v2_2; %Delta V at venus to put into parking orbit
magdV1_2 = norm(deltav1_2);
magdV2_2 = norm(deltav2_2);

disp(' ')
disp("--------Part 2--------")
disp("Magnitude of delta V1 (earth) : " + magdV1_2 + "km/s")
disp("Magnitude of delta V2 (venus capture): " + magdV2_2 + "km/s")
disp("Velcity upon entering Venus SOI: " + magv2inf2 + "km/s")


%------Venus to Mars-------
%
%

%posV and velV are still valid
[posM, velM] = planetEphemeris(juliandate(2005, 5, 14, 13, 23, 33), 'Sun', 'Mars', '405', 'km');
[posMPLOTTING, velMPLOTTING] = planetEphemeris(juliandate(2004, 6, 5, 1, 52, 21), 'Sun', 'Mars', '405', 'km');

ToF3 = 174.925; %days to go from earth to venus
[v1D3, v2A3] = lambert(posV, posM, ToF3, 0, 132712440018);


v1inf3 = v1D3 - velV; %v infinity departing Venus
v2inf3 = v2A3 - velM; % v infinity arriving at mars
magv1inf3 = norm(v1inf3); % magnitude of v1 infinity leaving venus
magv2inf3 = norm(v2inf3); % magnitude of v2 infinity arriving at mars

velop = sqrt(magv1inf3^2 + 2*324859/(4*rV)); %Velcoity of spacecraft at perigee on hyperbola at departure from venus
velop2 = sqrt(magv2inf3^2 + 2*42828/(4*rM)); %Velocity of spacecraft at pergiee on hyperbola at arrival at mars

vcircleMars = sqrt(42828/(4*rM));

deltav1_3 = velop - v2_2; %dvel = velocity on hyperbola at perigee - velocity on circular orbit.  Delta V at venus departure
deltav2_3 = velop2 - vcircleMars; %Delta V at mars to put into capture orbit
magdV1_3 = norm(deltav1_3);
magdV2_3 = norm(deltav2_3);


disp('--Venus to Mars--')
disp("Magnitude of delta V3 (venus departure) : " + magdV1_3 + "km/s")
disp("Magnitude of delta V4 (mars capture): " + magdV2_3 + "km/s")
disp("Velcity upon entering Mars SOI: " + magv2inf3 + "km/s")
z = v2A2 - v1D3;
fprintf('The vector is: [');
fprintf('%g ', z);
fprintf(']\n');

%-----PLOTTING-----
X = [posE,v1D2]; %initial parameters for spacecraft
time2 = 14467260;
timespan = [0,time2]; %Start and end times for ode45 to use
ode2=@drdr2;
options = odeset('AbsTol', 1e-9, 'RelTol', 1e-9); %Found to help keep tolerances from drifting values

%PLOTS TRANFER ORBIT EARTH TO VENUS
[t,y]= ode45(ode2,timespan,X, options); %Runs the ODE calculation numerically using a reference to the drdr function, time span, and initial conditions

xr = y(:,1); %Grabs x values of r
yr = y(:,2); %Grabs y values of r
zr = y(:,3); %Grabs z values of r

plot3(xr,yr,zr, 'linewidth', 2.5, 'Color', 'Black'); %Plots all values of r on 3d graph
hold on;

%PLOTS EARTH ORBIT
X2 = [posE, velE];

[t2,y2] = ode45(ode2, timespan, X2, options); %EARTH ORBIT DOESN'T PLOT ALL THE WAY, INTENTIONAL

xr2 = y2(:,1); %Grabs x values of r
yr2 = y2(:,2); %Grabs y values of r
zr2 = y2(:,3); %Grabs z values of r
plot3(xr2, yr2, zr2, 'linewidth', 2.5, 'Color', 'Blue')

%PLOTS SUN
plot3(0,0,0, '*', 'MarkerSize', 18, 'LineWidth', 2, 'Color', '#ffa500'); %Plots the sun

%PLOTS VENUS ORBIT
X3 = [posVPLOTTING, velVPLOTTING];
[t3,y3] = ode45(ode2, timespan, X3, options);

xr3 = y3(:,1); %Grabs x values of r
yr3 = y3(:,2); %Grabs y values of r
zr3 = y3(:,3); %Grabs z values of r
plot3(xr3, yr3, zr3, 'linewidth',2.5, 'Color', '#FFC649')



%PLOTS VENUS TO MARS TRANSFER
timed = [0,15113520]; %Time from venus to mars in seconds
X4 = [posV,v1D3]; %Initilal position and velocity of the satelite leaving venus
[t4,y4] = ode45(ode2, timed, X4, options);

xr4 = y4(:,1); %Grabs x values of r
yr4 = y4(:,2); %Grabs y values of r
zr4 = y4(:,3); %Grabs z values of r
plot3(xr4, yr4, zr4, 'linewidth',2.5, 'Color', 'Black')

%PLOTS MARS ORBIT
timer = [0, 29580768];
X5 = [posMPLOTTING, velMPLOTTING];
[t5,y5] = ode45(ode2, timer, X5, options);

xr5 = y5(:,1); %Grabs x values of r
yr5 = y5(:,2); %Grabs y values of r
zr5 = y5(:,3); %Grabs z values of r
plot3(xr5, yr5, zr5, 'linewidth',2.5, 'Color', 'Red')



view(0,89.9)
axis square
xlim([-3 * 10^8, 3 * 10^8])
ylim([-3 * 10^8, 3 * 10^8])







%--------------PART 3----------------

%Angle between incoming v_inf and outgoing v_inf = del 
del = acosd((dot(v2inf2, v1inf3))/(norm(v2inf2) * norm(v1inf3)));

%Solve for e
e = 1/sind(del/2);

%Solve for radius of perigee
rp = (e * 324859 - 324859) / (norm(v1inf3))^2;


dVv = v1inf3 - v2inf2;


%%------------FINAL TALLIES------------
disp(" ")
disp("---------Final-------")
zz = magdV1 + magdV2; %Part 1 total dv
disp("Part 1 total dV: " + zz + " km/s"); %total dv



z = v2A2 - v1D3; %DELTA V AT VENUS 
disp("Part 2 dV at Venus: " + norm(z) + " km/s")
total = norm(z) + magdV2_3 + magdV1_2; %TOtal part two dv
disp("Part 2 dv: " + total + " km/s")



disp("Part 3 dV from flyby: " + norm(dVv) + " km/s");
total2 = magdV1_2 + magdV2_3;
disp("Part 3 total mission dV: " + total2 + " km/s")

disp("Perigee radius of Venus flyby: " + rp + " km")
figure();

h = animatedline("Color", 'Black');
h2 = animatedline('Color','Blue');
h3 = animatedline('Color','#FFC649');
h4 = animatedline('Color','Black');
h5 = animatedline('Color','Red');
view(3);

timer = [0, 29580768];
X = [posE,v1D2]; %initial parameters for spacecraft
time2 = 14467260;
timespan = [0,timer]; %Start and end times for ode45 to use
ode2=@drdr2;
options = odeset('AbsTol', 1e-9, 'RelTol', 1e-9); %Found to help keep tolerances from drifting values

%PLOTS TRANFER ORBIT EARTH TO VENUS
[t,y]= ode45(ode2,[0,time2],X, options); %Runs the ODE calculation numerically using a reference to the drdr function, time span, and initial conditions

xr = y(:,1); %Grabs x values of r
yr = y(:,2); %Grabs y values of r
zr = y(:,3); %Grabs z values of r


hold on;

%PLOTS EARTH ORBIT
X2 = [posE, velE];

[t2,y2] = ode45(ode2, timer, X2, options); %EARTH ORBIT DOESN'T PLOT ALL THE WAY, INTENTIONAL

xr2 = y2(:,1); %Grabs x values of r
yr2 = y2(:,2); %Grabs y values of r
zr2 = y2(:,3); %Grabs z values of r


%PLOTS SUN


%PLOTS VENUS ORBIT
X3 = [posVPLOTTING, velVPLOTTING];
[t3,y3] = ode45(ode2, timer, X3, options);

xr3 = y3(:,1); %Grabs x values of r
yr3 = y3(:,2); %Grabs y values of r
zr3 = y3(:,3); %Grabs z values of r




%PLOTS VENUS TO MARS TRANSFER
timed = [0,15113520]; %Time from venus to mars in seconds
X4 = [posV,v1D3]; %Initilal position and velocity of the satelite leaving venus
[t4,y4] = ode45(ode2, timed, X4, options);

xr4 = y4(:,1); %Grabs x values of r
yr4 = y4(:,2); %Grabs y values of r
zr4 = y4(:,3); %Grabs z values of r


%PLOTS MARS ORBIT
timer = [0, 29580768];
X5 = [posMPLOTTING, velMPLOTTING];
[t5,y5] = ode45(ode2, timer, X5, options);

xr5 = y5(:,1); %Grabs x values of r
yr5 = y5(:,2); %Grabs y values of r
zr5 = y5(:,3); %Grabs z values of r

xrfinal = [xr; xr4];
yrfinal = [yr; yr4];
zrfinal = [zr; zr4];

tfinal = [t;t4+max(t)];
view(0,89.9)
j = 2.5*10^8
xlim([-j,j]);
ylim([-j,j]);

%index1 = 0;
index2 = 1;
index3 = 1;
%index4 = 0;
index5 = 1;
indexfinal = 1;
n = 500;
for k = 1:n;

    ixr2 = round(k/n * length(xr2));
    ixr3 = round(k/n * .8*length(xr3));
    ixr5 = round(k/n * length(xr5));
    ixrfinal = round(k/n * length(xrfinal));

    if(ixr5<1)
        ixr5 = 1;
    end
    addpoints(h2,xr2(ixr2),yr2(ixr2), zr2(ixr2));
    addpoints(h3,xr3(ixr3),yr3(ixr3), zr3(ixr3));
    addpoints(h5,xr5(ixr5),yr5(ixr5), zr5(ixr5));
    addpoints(h4, xrfinal(ixrfinal), yrfinal(ixrfinal), zrfinal(ixrfinal))
    drawnow 
end

% figur




