function P = proRobGetCoords(robot)
%PROROBGETCOORDS - get coordinates of the robot
%
% Args:
% robot - robot handle
%
% Return coorinates in mm and angles (roll, pitch, yaw) in degrees

	%just check if it's not virtual robot :P
	if (~isa(robot, 'mmRobotDef') || ~strcmp(robot.comtype, 'Serial'))
		error('Error reading the location');
	end

	J = mmGetPos2(robot, 'J');	
	
	[m n] = size(J);
	if (m ~= 6 || n ~= 1)
		error('Could not read joint coordinates from the robot.');
	end

	%%% ANGULAR MISALIGNMENT %%%
	AlphaOffRad = [ 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000 ];
	%%% MISALIGNMENT %%%
	aMis = [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ];
	dMis = [ 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 ];
	%%% END EFFECTOR MISALIGNMENT - TOOL %%%
	toolMis = [ 0.00; 0.00; 0.00 ];

	DH = robot.denavitHartenberg;
	DH(1,:) = DH(1,:) + AlphaOffRad;
	DH(2,:) = DH(2,:) + aMis;
	DH(4,:) = DH(4,:) + dMis;

	%compute DKT
	J = J*pi/180; %deg to rad
	T = robot.base;

	for i=1:size(DH, 2)
		T = T * denavitHartenberg(DH(:,i) + J(i) * robot.denavitHartenbergParameters(:,i));
	end

	T = T * [1 0 0 toolMis(1); 0 1 0 toolMis(2); 0 0 1 toolMis(3); 0 0 0 1];
	X = T * [0; 0; 0; 1];
	[A B C] = mtxToAngles(T);

	P = [X(1); X(2); X(3); A; B; C];


end

