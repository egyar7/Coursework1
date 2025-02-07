%Adil Rafique 
%egyar7@nottingham.ac.uk

%% Q1 - VARIABLES [8 MARKS]
clear

%Defining x:
x = 1.5;

% a)
a = sin(2*x);
% b)
b = sin(x)^2;
% c)
c = x^5 + 3*x^4 + x^3 + (x^2)/2 + x;
% d)
d = exp(sin(x)*cos(x));
% e)
e = (x + 1)*(x-1)*((x^2) + 1)*(x^3);
% f)
f = a + 2*b;
% g)
g = (f/(c + e))^2;

%% Q2 - VECTORS [12 MARKS]
clear

%Defining x:
x = linspace(-2*pi,2*pi,100);

% a)
a = sin(x);
% b)
b = cos(x).^2;
% c)
c = tan(x);
% d)
d = sinh(x);
% e)
average = mean(b); 
standardDeviation = std(b);
% f) 
% The mean is not exactly 5 because only 100 values were taken and
% approximated. Using 100 points to approximate the value of the continuous 
% function means that it cannot capture the whole function accurately. 
% If the number of values was increased, the mean would be closer to the
% true answer as it would be a better approximation.

% g)
%Defining the waves time vector
timeVector = linspace(0,2*pi,100);

%Constants from Question
iPeak = 0.1;
resistance = 1000;

%Calculating the current at each point
currentValues = iPeak * sin(timeVector);

%Calculating the power at each point
powerValues = currentValues.^2 * resistance;

%Finding average power
averagePower = mean(powerValues);

fprintf('The average power is %.4f Watts\n',averagePower)

%% Q3 - SOLVING EQUATIONS [6 MARKS]
clear

% a)
%polynomial = 24*x^2 - 24*x - 480;

%Defining the coefficients
a = 24;
b = -24;
c = -480;

%Substitution into quadratic formula to find roots.
x1 = (-b+sqrt(b^2-4*a*c))/(2*a);
x2 = (-b-sqrt(b^2-4*a*c))/(2*a);

fprintf('The roots of the polynomial occur when x = %.2f and x = %.2f\n',x1,x2);

% bi)

%Defining constants
epsilon0 = 8.85418782e-12;
mu0 = (4 * pi) * 1e-7;

%Calculation using formula
speedOfLight = sqrt(1/(epsilon0*mu0));

fprintf('The speed of light is %f m/s.\n',speedOfLight)

% bii)

%Defining constants
lowerWaveLength = 0.6e-6;
upperWaveLength = 28e-6;

%Calculating the frequency
lowerWaveLengthToFrequency = speedOfLight/lowerWaveLength;
upperWaveLengthToFrequency = speedOfLight/upperWaveLength;

fprintf('The range of electromagnetic frequencies are %.2f Hz and %.2fHz.\n',lowerWaveLengthToFrequency,upperWaveLengthToFrequency)

%% Q4 - EXAMPLES [12 MARKS]
clear
