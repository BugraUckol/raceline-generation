function C = CE2B(E)

PHI = E(1);
THE = E(2);
PSI = E(3);

% Transformation between the Earth frame and first intermediate frame
% Rotating about z-axis for PSI radians
CE2F = [ cos(PSI), sin(PSI), 0;
        -sin(PSI), cos(PSI), 0;
                0,        0, 1];

% Transformation between the first and second intermatediate frames
% Rotating about y-axis for THE (theta) radians
CF2S = [ cos(THE), 0, -sin(THE);
                0, 1,         0;
         sin(THE), 0,  cos(THE)];

% Transformation between the second intermediate frame and body frame
% Rotating about x-axis for PHI radians
CS2B = [ 1,         0,        0;
         0,  cos(PHI), sin(PHI);
         0, -sin(PHI), cos(PHI)];

% Combination
C = CS2B * CF2S * CE2F;

end