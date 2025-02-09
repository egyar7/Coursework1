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

% a) 
% Clear removes the variables that are stored in the memory.

% b)
% Who shows a list of all the variables that are in the current workspace.

% c) 
% A semicolon after a line prevents the code in that line from being
% displayed.

% d)
% Sound plays the sound relating to a vector with a specific sampling rate.

% e)
% Roots finds the roots of a polynomial when given the coefficients.

% f)
% Abs finds the absolute value of a number or an array. 

%% Q5 - PROGRAM FLOW [12 MARKS]
clear,clc,clearvars

% a) Flowchart done in different program

% b and c) 

%Variables

initialHeight = 10000;
time = 0;
velocity = 200;
timeInterval = 2;

height = initialHeight;

parachuteDeployment = false;

%Height Verification 

while height >= 0

    %Loop to activate parachute

    if height <= 2700 && ~parachuteDeployment
        disp('Parachute Deployed');
        disp(['Time elapsed before parachute deployment: ', num2str(time), ' seconds']);
        parachuteDeployment = true;
    end
    
    %Loop to activate transponder beacon

    if height <= 0
        disp('Transponder Beacon Activated');
        break;
    end
    
    %Decreasing the height per time interval

    height = height - velocity * timeInterval;

    %Increasing the time interval

    time = time + timeInterval;

end

%% Q6 FORMAT AND PRINT TEXT TO SCREEN AND LOOPS [13 MARKS]
clear,clc,clearvars

%Defining arrays
time = [1300 , 1600 , 1900];
temperature = [19 , 20 , 18];
humidity = [55 , 49 , 59];
uvLevel = [4 , 2 , 1];

%Displaying the date and location
disp('Data logging initiatied - 27/7/2023')
disp('Location - Nottingam')
disp(' ')

%Loop that displays the weather information
%(Pulls a value from the next column in each array after every run)
%(until each value in the array has been used.)
n = 1;
while n <= length(time)
    statement = sprintf('Time\t\t%d\nTemperature\t%d C\nHumidity\t%d%%\nUV level\t%d\n', ...
        time(n),temperature(n),humidity(n),uvLevel(n));
    disp(statement);
    n = n+1;
end

disp('Data logging terminated')

%The limitation is that the data has to be input into the code for it to be
%displayed. This means that it cannot get realtime weather information
%unless it is constantly input into the code. To improve the program, the
%weather data should be input automatically from the sensors into the
%code by reading data provided by the sensors.

%% Q7 - FOR LOOPS AND DISPLAYING DATA [16 MARKS]
clear,clc,clearvars

%Defining number of terms
numberOfTerms = 50;

%Preallocating array for the fibonacci sequence
fibonacci = zeros(1,numberOfTerms);

%Setting the first and second term of the fibonacci sequence
fibonacci(1) = 0;
fibonacci(2) = 1;

%Preallocating array for the golden ratio
goldenRatio = zeros(1,numberOfTerms);

%For loop to generate fibonacci sequence and golden ratio
for n = 3:numberOfTerms
    fibonacci(n) = fibonacci(n-1) + fibonacci(n-2);
    goldenRatio(n) = fibonacci(n)/fibonacci(n-1);
end

%Statement and headings for the table of values.
fprintf('Fibonacci sequence and approximations to the golden ratio:\n');
fprintf('%-5s %-20s %-20s\n', 'Index', 'Fibonacci Number', 'Golden Ratio');

%For loop to display the fibonacci sequence and corresponding golden ratio.
for n = 1:50
    if n > 1
        fprintf('%-5d %-20d %-20.5f\n', n, fibonacci(n), goldenRatio(n));
    else
        
        %First golden ratio cannot be calculated so 'N/A' must display.
        fprintf('%-5d %-20d %-20s\n', n, fibonacci(n), 'N/A');
    end
    
    %Convergence function to end loop when golden ratios are within 0.1%.
    if n > 1 && abs((goldenRatio(n)-goldenRatio(n-1))/goldenRatio(n-1)) <= 0.001
        disp('Golden ratio converged to within 0.1%');
        break;
    end
end

%Part d stops the loop when the golden ratio converges. This removes
%unnecessary outputs and makes the code more efficient. It also reduces the
%amount of memory used and increases the performance/speed at which the
%program is completed.

%% Q8 - USING THE SWITCH STATEMENT [10 MARKS]
clear

