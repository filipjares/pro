function r = proServoOff(robot)
%PROSERVOOFF - turns robot servo off
%
% Args:
% robot - robot handle
%
% Returns 0 for success.

	%turn off servo
	r = mmSrvOff(robot);
	
end

