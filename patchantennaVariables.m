clear; close all; clc
prompt = "What is the desired frequency (GHz) =  ";
f = input(prompt);
prompt = "What is the substrate's dielectric constant =  ";
Er = input(prompt);
prompt = "What is the substrate's height (mm) =  ";
h = input(prompt);
prompt = "What is the characteristic impedence (ohm) =  ";
Zc = input(prompt);
h=h/1000;
c=299792458;
%%Calculating the Width
sqrt_i=2/(Er+1);
W=((c)/(2*f*10^9))*(sqrt(sqrt_i));
fprintf("Width of patch = %.2f mm ",W*1000);
%%effective dielectiric constant
sqrt_i=1+12*h/W;
Eeff=((Er+1)/2)+((Er-1)/2)*(1/sqrt(sqrt_i));
%%fprintf("\nEeff = %.2f",Eeff);
%%Extension Length
dl=(h*0.412)*((Eeff+0.3)/(Eeff-0.258))*((W/h+0.284)/(W/h+0.8));
%%Actual Length
L=(c)/(2*f*10^(9)*sqrt(Eeff))-(2*dl);
fprintf("\nLength of patch = %.2f mm ",L*1000);
%%Gap length
g=46*c/(2*Eeff*f*10^9);
fprintf("\nLength of gap = %.2f mm ",g);
% Width of the Feed Line (W0)
h=h*1000;
A = (Zc/60)*sqrt((Er + 1)/2) + ((Er - 1)/(Er + 1))*(0.23 + 0.11/Er);
B = ((377*pi)/(2*Zc*sqrt(Er)));
W_A = (8*exp(A))/(exp(2*A)-2);
W_B = (2/pi)*(B-1-log(2*B-1) + ((Er - 1)/(2*Er))*(log(B-1) + 0.39 - 0.61/Er) );
if W_A < 2
    W0 = W_A*h; % in mm
elseif W_B > 2
    W0 = W_B*h; % in mm
end
fprintf("\nWidth of the feedline = %.2f mm ",W0);


% Characteristic Impedence of the feed line (Z0) 
% In the calculation we have taken Z0 as 50 (Z0=Zc). But also we can
% calculate the real Z0 as given below.
if (W0/h) <= 1
    Z0 = (60/sqrt(Eeff))*log((8*h*10^(-3)/W0) + (W0/(4*h*10^(-3))));
elseif (W0/h) > 1
    Z0 = (120*pi/sqrt(Eeff))/((W0/h) + 1.393 + 0.667*log((W0/h) + 1.444));
end



% Antenna Input Impedence (Rin)
k0 = (2*pi*f*10^9/c);
x = k0*W ; 
Si_x = sinint(x) ; % Sine integeral
I_1 = -2 + cos(x) + x*Si_x + sin(x)/x ;
G_1 = I_1/(120*pi.^2);
fun = @(theta) (((sin(k0*W*cos(theta)/2))/(cos(theta)))^2).*((sin(theta)).^3).*besselj(1,k0*L*sin(theta)) ; 
G_12 = (1/(120*pi.^2))*integral(fun,0,pi);
Rin_odd = 1/(2*(G_1 + G_12)) ; % in Ohms (for odd modes)
Rin_even = 1/(2*(G_1 - G_12)) ; % in Ohms (for even modes)
% Position of the inset-feed point (y0)
y0 = acos(abs(sqrt(Zc/Rin_odd)))*L/pi ; % in meters

fprintf('\nDepth of Inset-feed (y0) = %.3f mm\n',y0*1000);
L=L*1000;
W=W*1000;
y0=y0*1000;

fprintf("           Dimensions          \n//////////////////////////////////////////////\n\n\n");
fprintf("Ground\n\n");
fprintf("Position X = %.2f    Y = %.2f   Z = %.3f  (mm) \n",0,0,0);
fprintf("Size     X = %.2f    Y = %.2f   Z = %.3f  (mm) \n\n",2*L,3*W,0.035);
fprintf("Substrate\n\n");
fprintf("Position X = %.2f    Y = %.2f   Z = %.3f  (mm) \n",0,0,0.035);
fprintf("Size     X = %.2f    Y = %.2f   Z = %.3f  (mm) \n\n",2*L,3*W,h);
fprintf("Patch\n\n");
fprintf("Position X = %.2f    Y = %.2f   Z = %.3f  (mm) \n",(L/2),W,h+0.035);
fprintf("Size     X = %.2f    Y = %.2f   Z = %.3f  (mm) \n\n",L,W,0.035);
fprintf("Gap\n\n");
fprintf("Position X = %.2f    Y = %.2f   Z = %.3f  (mm) \n",(L-y0)+L/2,W+(W/2)-g-(W0/2),h+0.035);
fprintf("Size     X = %.2f    Y = %.2f   Z = %.3f  (mm) \n\n",y0,W0+2*g,0.035);
fprintf("Transmission Line\n\n");
fprintf("Position X = %.2f    Y = %.2f   Z = %.3f  (mm) \n",(L-y0)+L/2,W+(W/2)-(W0/2),h+0.035);
fprintf("Size     X = %.3f    Y = %.2f   Z = %.3f  (mm) \n\n",y0+L/2,W0,0.035);
fprintf("Source\n\n");
fprintf("Position X = %.3f    Y = %.2f   Z = %.3f  (mm) \n",2*L,W+(W/2)-(W0/2),0.035);
fprintf("Size                 Y = %.2f   Z = %.2f  (mm) \n\n",W0,h+0.035);



