function robot = proRobOpen()
%PROROBOPEN - opens robot handle
%
% Returns robot handle
% jCoords - joint coordinates


	%open robot with procedural serial line
	robot = mmOpen('RV6S');
	%props.instantCom = 1;
	%props.instantMov = 1;
	%robot = mmOpenVirt('RV6S', 'comprops', props);

	%control on
	r = mmCntlOn(robot);
		
	if (r)
		
		%quit if epic fail
		mmClose(robot);
		error('Failed opening the robot.');
		
	end
	
end

