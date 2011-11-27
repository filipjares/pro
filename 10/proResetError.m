function r = proResetError( robot )
%PRORESETERROR - resets robot error
%
% Args:
% robot - robot handle
%
% Returns 0 for success.

	%purge error
	r = mmPurge(robot);

end

