function C = CB2W(ALPHA, BETA)

% ALPHA HAS NEGATIVE SENSE BY DEFINITION!
% Transformation between the first and second intermatediate frames
% Rotating about y-axis for ALPHA radians
CF2S = [  cos(ALPHA), 0,  sin(ALPHA);
                   0, 1,           0;
         -sin(ALPHA), 0,  cos(ALPHA)];

% Transformation between the second intermediate frame and body frame
% Rotating about z-axis for BETA radians
CS2B = [ cos(BETA), sin(BETA), 0;
        -sin(BETA), cos(BETA), 0;
                 0,         0, 1];

% Combination
C = CS2B * CF2S;

end