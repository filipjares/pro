function r = proRobMove( robot, jCoords )
%PROROBMOVE - moves robot to joint coordinates (requires servo on)
%
% Args:
% robot - robot handle
% jCoords - joint coordinates
%
% Returns 0 for success.

	%check coords
	[m n] = size(jCoords);
	
	if (m ~= 6 || n ~= 1)
		error('Invalid J Coodinates');
	end

	%just check if it's not virtual robot :P
	if (~isa(robot, 'mmRobotDef') || ~strcmp(robot.comtype, 'Serial'))
		error('Error moving the robot');
	end
	
	%%% ANGULAR OFFSETS %%%
	ThetaOffDeg = [ 13.9451; 11.1378; 13.2438; 7.5678; 15.5907; -278.8540 ];
	jCoords = jCoords + ThetaOffDeg;
	
	r = mmMovSafe(robot, 'J', jCoords);

end

