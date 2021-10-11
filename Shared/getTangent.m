function t = getTangent(jacobian)
%getTangent
% it is assumed that jacobian has full rank
    A     = jacobian;
    dim   = size(A,2)-size(A,1); % dimension of null space
    if dim > 0
        [Q,~] = qr(A');
        z     = Q(:,(end+1-dim):end); % tangent space
        if det([A;z'])<0 && dim == 1
            t     = -z;
        else
            t     = z;
        end
    elseif dim == 0
        t = [];
    else
        error('jacobian is of wrong dimension')
    end
end