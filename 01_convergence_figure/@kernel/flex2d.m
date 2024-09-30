function obj = flex2d(type, zk, coefs)
%KERNEL.FLEX2D   Construct the Flexural wave kernels for evaluations.

%
% See CHNK.FLEX2D.KERN.

if ( nargin < 1 )
    error('Missing kernel type.');
end

if ( nargin < 2 )
    error('Missing wavenumber.');
end

obj = kernel();
obj.name = 'flexural';
obj.params.zk = zk;
obj.opdims = [1 1];

switch lower(type)

    case{'clamped-plate-first', 'first-kernel'}
        obj.type = 'clamped-plate-first';
        obj.eval = @(s,t) chnk.flex2d.kern(zk, s, t, 'first kernel');
      
        obj.sing = 'log';
    case{'clamped-plate-second', 'second-kernel'}
        obj.type = 'clamped-plate-second';
        obj.eval = @(s,t) chnk.flex2d.kern(zk, s, t, 'second kernel');
      
        obj.sing = 'log';

    case{'free-plate-eval-first', 'free plate first kernel'}
        obj.type = 'free-plate-eval-first';
        obj.eval = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first');
     
        obj.sing = 'log';

    case{'free-plate-eval-first hilbert', 'free plate first kernel hilbert'}
        obj.type = 'free-plate-eval-first-hilbert';
        obj.eval = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval first hilbert', coefs);
     
        obj.sing = 'log';

    case{'free-plate-eval-second', 'free plate second kernel'}
        obj.type = 'free-plate-eval-second';
        obj.eval = @(s,t) chnk.flex2d.kern(zk, s, t, 'free plate eval second');
      
        obj.sing = 'log';     

    otherwise
        error('Unknown flexural wave kernel type ''%s''.', type);

end

end




