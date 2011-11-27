function r = proServoOn(robot)
%PROSERVOON - turns robot servo on
%
% Args:
% robot - robot handle
%
% Returns 0 for success.

	%turn on servo
	r = mmSrvOn(robot);
	
end