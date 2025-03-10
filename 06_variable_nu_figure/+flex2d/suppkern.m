function submat= suppkern(zk,srcinfo,targinfo,type,varargin)
% suppkern Modified biharmonic layer potential kernels for the supported
% plate
% 
% Syntax: submat = suppkern(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
%
% Input:
%   zk - complex number, Helmholtz wave number
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   targinfo - description of targets in ptinfo struct format,
%                if info not relevant (d/d2) it doesn't need to
%                be provided. sprime requires tangent info in
%                targinfo.d
%   type - string, determines kernel type
%                type == 'd', double layer kernel D
%                type == 's', single layer kernel S
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'dprime', normal derivative of double layer D'
%                type == 'c', combined layer kernel coef(1) D + coef(2) S
%                type == 'stau', tangential derivative of single layer
%                type == 'all', returns all four layer potentials, 
%                       [coef(1,1)*D coef(1,2)*S; coef(2,1)*D' coef(2,2)*S']
%                type == 'c2trans' returns the combined field, and the 
%                          normal derivative of the combined field
%                        [coef(1)*D + coef(2)*S; coef(1)*D' + coef(2)*S']
%                type == 'trans_rep' returns the potential corresponding
%                           to the transmission representation
%                        [coef(1)*D coef(2)*S]
%                type == 'trans_rep_prime' returns the normal derivative
%                          corresponding to the transmission representation
%                        [coef(1)*D' coef(2)*S']
%                type == 'trans_rep_grad' returns the gradient corresponding
%                         to the transmission representation
%                        [coef(1)*d_x D coef(2)*d_x S;
%                         coef(1)*d_y D coef(2)*d_y S]
%
%   varargin{1} - coef: length 2 array in the combined layer 
%                 formula, 2x2 matrix for all kernels
%                 otherwise does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also FLEX2D.KERN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

% supported plate kernel for modified biharmonic problem
if strcmpi(type, 'supported plate')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   srcd2 = srcinfo.d2;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};
   nu = coefs(1);

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   dx1 = repmat(srctang(1,:),nt,1);
   dy1 = repmat(srctang(2,:),nt,1);
   
   d2x1 = repmat(srcd2(1,:),nt,1);
   d2y1 = repmat(srcd2(2,:),nt,1);
   
   denom = sqrt(dx1.^2+dy1.^2).^3;
   numer = dx1.*d2y1-d2x1.*dy1;

   kappa = numer./denom;   

   xs = repmat(src(1,:),nt,1);
   ys = repmat(src(2,:),nt,1);

   % Ellipse:  

   c = 1;

   t = atan2(c*ys,xs);

   % xs = c*cos(t); x coordinate of parametrization
   % ys = sin(t); y coordinate of parametrization

   x1 = -c*sin(t);    
   y1 = cos(t);

   x2 = -c*cos(t);
   y2 = -sin(t);

   x3 = c*sin(t);
   y3 = -cos(t);

   % % Starfish :  
   % 
   % narms = 3;
   % amp = 0.3;
   % phi = pi/3;
   % 
   % t = atan2(ys,xs);
   % 
   % ct = cos(t);
   % st = sin(t);
   % cnt = cos(narms*(t + phi));
   % snt = sin(narms*(t + phi));
   % 
   % x1 = -(1+amp*cnt).*st-narms*amp*snt.*ct;
   % y1 = (1+amp*cnt).*ct-narms*amp*snt.*st;
   % 
   % x2 = -y1-narms*amp*(narms*cnt.*ct-snt.*st);
   % y2 = x1-narms*amp*(narms*cnt.*st+snt.*ct);
   % 
   % x3 = -y2 + amp*narms^3*snt.*ct + 2*amp*narms^2*cnt.*st + amp*narms*snt.*ct;
   % y3 = x2 + amp*narms^3*snt.*st - 2*amp*narms^2*cnt.*ct + amp*narms*snt.*st;
    
   N = x1.*y2 - y1.*x2;
   dNdt = x1.*y3 - y1.*x3;
   D = (x1.^2 + y1.^2).^( 3/2);
   dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
   kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2);

   a1 = 2-nu;
   a2 = (-1+nu)*(7+nu)/(3 - nu);
   a3 = (1-nu)*(3+nu)/(1+nu);
   
   [~, grad, hess, third, forth, ~] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part

   [~, ~, ~, ~, ~, fifth] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part

   K11 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

   K21 = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % + ... % G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % nu*G_{taux taux tauy}

   K12 = -1/(2*zk^2).*(grad(:, :, 1).*nx + grad(:, :, 2).*ny);
           
   K22 = -1/(2*zk^2).*(third(:,:,1).*(nx.*nxtarg.^2) + third(:,:,2).*(2*nx.*nxtarg.*nytarg + ny.*nxtarg.^2) + ...
         + third(:,:,3).*(nx.*nytarg.^2 + 2*ny.*nxtarg.*nytarg) + third(:,:,4).*(ny.*nytarg.^2))+...
         -nu./(2*zk^2).*(third(:,:,1).*(nx.*tauxtarg.^2) + third(:,:,2).*(2*nx.*tauxtarg.*tauytarg + ny.*tauxtarg.^2) + ...
         + third(:,:,3).*(nx.*tauytarg.^2 + 2*ny.*tauxtarg.*tauytarg) + third(:,:,4).*(ny.*tauytarg.^2)); 

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21; 
  submat(2:2:end,2:2:end) = K22; 

end

% supported plate kernel for modified biharmonic problem
if strcmpi(type, 'supported plate ellipse')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   srcd2 = srcinfo.d2;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};
   nu = coefs(1);

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   dx1 = repmat(srctang(1,:),nt,1);
   dy1 = repmat(srctang(2,:),nt,1);
   
   d2x1 = repmat(srcd2(1,:),nt,1);
   d2y1 = repmat(srcd2(2,:),nt,1);
   
   denom = sqrt(dx1.^2+dy1.^2).^3;
   numer = dx1.*d2y1-d2x1.*dy1;

   kappa = numer./denom;   

   xs = repmat(src(1,:),nt,1);
   ys = repmat(src(2,:),nt,1);

   % Ellipse:  

   c = 2;

   t = atan2(c*ys,xs);

   x1 = -c*sin(t);    
   y1 = cos(t);

   x2 = -c*cos(t);
   y2 = -sin(t);

   x3 = c*sin(t);
   y3 = -cos(t);

   N = x1.*y2 - y1.*x2;
   dNdt = x1.*y3 - y1.*x3;
   D = (x1.^2 + y1.^2).^( 3/2);
   dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
   kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2); % kappa prime

   a1 = 2-nu;
   a2 = (-1+nu)*(7+nu)/(3 - nu);
   a3 = (1-nu)*(3+nu)/(1+nu);
   
   [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);      

   K11 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

   K12 = -1/(2*zk^2).*(grad(:, :, 1).*nx + grad(:, :, 2).*ny);
           
   K22 = -1/(2*zk^2).*(third(:,:,1).*(nx.*nxtarg.^2) + third(:,:,2).*(2*nx.*nxtarg.*nytarg + ny.*nxtarg.^2) + ...
         + third(:,:,3).*(nx.*nytarg.^2 + 2*ny.*nxtarg.*nytarg) + third(:,:,4).*(ny.*nytarg.^2))+...
         -nu./(2*zk^2).*(third(:,:,1).*(nx.*tauxtarg.^2) + third(:,:,2).*(2*nx.*tauxtarg.*tauytarg + ny.*tauxtarg.^2) + ...
         + third(:,:,3).*(nx.*tauytarg.^2 + 2*ny.*tauxtarg.*tauytarg) + third(:,:,4).*(ny.*tauytarg.^2)); 


   [~, ~, ~, third, forth, fifth] = flex2d.hkdiffgreen(zk, src, targ, true);      

   K21 = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ... % K21 with the biharmonic Green's function subtracted off (biharmonic part handled in a separate function)
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21; 
  submat(2:2:end,2:2:end) = K22; 

end


% supported plate kernel for modified biharmonic problem
if strcmpi(type, 'supported plate droplet')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    targnorm = targinfo.n;
    targtang = targinfo.d;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    
    zkimag = (1i)*zk;
    
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    rx = targ(1,:).' - src(1,:);
    ry = targ(2,:).' - src(2,:);
    r2 = rx.^2 + ry.^2;
    
    dx1 = repmat((targtang(1,:)).',1,ns);
    dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    tauxtarg = dx1./ds1;
    tauytarg = dy1./ds1;
    
    dx1 = repmat(srctang(1,:),nt,1);
    dy1 = repmat(srctang(2,:),nt,1);
    
    d2x1 = repmat(srcd2(1,:),nt,1);
    d2y1 = repmat(srcd2(2,:),nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    kappa = numer./denom;   
    
    xs = repmat(src(1,:),nt,1);
    ys = repmat(src(2,:),nt,1);
    
    % droplet: 
    
    a = 2;
    b = 1;
    c = -0.4;
    
    t = atan2(a*ys - c*xs.^2/a, b*xs);
    
    % x = a*cos(t);
    % y = b*sin(t) + c*cos(t).^2;
    
    x1 = -a*sin(t);
    y1 = b*cos(t) - 2*c*cos(t).*sin(t);
    
    x2 = -a*cos(t);
    y2 = -b*sin(t) + 2*c*sin(t).^2 - 2*c*cos(t).^2;
    
    x3 = a*sin(t);
    y3 = -b*cos(t) + 8*c*sin(t).*cos(t);
    
    N = x1.*y2 - y1.*x2;
    dNdt = x1.*y3 - y1.*x3;
    D = (x1.^2 + y1.^2).^( 3/2);
    dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
    kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2); % kappa prime
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);      
    
    K11 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);
    
    K12 = -1/(2*zk^2).*(grad(:, :, 1).*nx + grad(:, :, 2).*ny);
           
    K22 = -1/(2*zk^2).*(third(:,:,1).*(nx.*nxtarg.^2) + third(:,:,2).*(2*nx.*nxtarg.*nytarg + ny.*nxtarg.^2) + ...
         + third(:,:,3).*(nx.*nytarg.^2 + 2*ny.*nxtarg.*nytarg) + third(:,:,4).*(ny.*nytarg.^2))+...
         -nu./(2*zk^2).*(third(:,:,1).*(nx.*tauxtarg.^2) + third(:,:,2).*(2*nx.*tauxtarg.*tauytarg + ny.*tauxtarg.^2) + ...
         + third(:,:,3).*(nx.*tauytarg.^2 + 2*ny.*tauxtarg.*tauytarg) + third(:,:,4).*(ny.*tauytarg.^2)); 
    
    
    [~, ~, ~, third, forth, fifth] = flex2d.hkdiffgreen(zk, src, targ, true);      
    
    K21 = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ... % K21 with the biharmonic Green's function subtracted off (biharmonic part handled in a separate function)
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}
    
    submat = zeros(2*nt,2*ns);
    
    submat(1:2:end,1:2:end) = K11;
    submat(1:2:end,2:2:end) = K12;
    
    submat(2:2:end,1:2:end) = K21; 
    submat(2:2:end,2:2:end) = K22; 

end


% supported plate kernel for modified biharmonic problem
if strcmpi(type, 'supported plate starfish 3arms 1')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    targnorm = targinfo.n;
    targtang = targinfo.d;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    
    zkimag = (1i)*zk;
    
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    rx = targ(1,:).' - src(1,:);
    ry = targ(2,:).' - src(2,:);
    r2 = rx.^2 + ry.^2;
    
    dx1 = repmat((targtang(1,:)).',1,ns);
    dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    tauxtarg = dx1./ds1;
    tauytarg = dy1./ds1;
    
    dx1 = repmat(srctang(1,:),nt,1);
    dy1 = repmat(srctang(2,:),nt,1);
    
    d2x1 = repmat(srcd2(1,:),nt,1);
    d2y1 = repmat(srcd2(2,:),nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    kappa = numer./denom;   
    
    xs = repmat(src(1,:),nt,1);
    ys = repmat(src(2,:),nt,1);
    
    % Starfish :  

    narms = 3;
    amp = 0.3;
    phi = 0;

    t = atan2(ys,xs);

    ct = cos(t);
    st = sin(t);
    cnt = cos(narms*(t + phi));
    snt = sin(narms*(t + phi));

    x1 = -(1+amp*cnt).*st-narms*amp*snt.*ct;
    y1 = (1+amp*cnt).*ct-narms*amp*snt.*st;

    x2 = -y1-narms*amp*(narms*cnt.*ct-snt.*st);
    y2 = x1-narms*amp*(narms*cnt.*st+snt.*ct);

    x3 = -y2 + amp*narms^3*snt.*ct + 2*amp*narms^2*cnt.*st + amp*narms*snt.*ct;
    y3 = x2 + amp*narms^3*snt.*st - 2*amp*narms^2*cnt.*ct + amp*narms*snt.*st;
    
    N = x1.*y2 - y1.*x2;
    dNdt = x1.*y3 - y1.*x3;
    D = (x1.^2 + y1.^2).^( 3/2);
    dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
    kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2); % kappa prime
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);      
    
    K11 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);
    
    K12 = -1/(2*zk^2).*(grad(:, :, 1).*nx + grad(:, :, 2).*ny);
           
    K22 = -1/(2*zk^2).*(third(:,:,1).*(nx.*nxtarg.^2) + third(:,:,2).*(2*nx.*nxtarg.*nytarg + ny.*nxtarg.^2) + ...
         + third(:,:,3).*(nx.*nytarg.^2 + 2*ny.*nxtarg.*nytarg) + third(:,:,4).*(ny.*nytarg.^2))+...
         -nu./(2*zk^2).*(third(:,:,1).*(nx.*tauxtarg.^2) + third(:,:,2).*(2*nx.*tauxtarg.*tauytarg + ny.*tauxtarg.^2) + ...
         + third(:,:,3).*(nx.*tauytarg.^2 + 2*ny.*tauxtarg.*tauytarg) + third(:,:,4).*(ny.*tauytarg.^2)); 
    
    
    [~, ~, ~, third, forth, fifth] = flex2d.hkdiffgreen(zk, src, targ, false);      
    
    K21 = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ... % K21 with the biharmonic Green's function subtracted off (biharmonic part handled in a separate function)
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}
    
    submat = zeros(2*nt,2*ns);
    
    submat(1:2:end,1:2:end) = K11;
    submat(1:2:end,2:2:end) = K12;
    
    submat(2:2:end,1:2:end) = K21; 
    submat(2:2:end,2:2:end) = K22; 

end


% supported plate kernel for modified biharmonic problem
if strcmpi(type, 'supported plate starfish 3arms 2')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    targnorm = targinfo.n;
    targtang = targinfo.d;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    
    zkimag = (1i)*zk;
    
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    rx = targ(1,:).' - src(1,:);
    ry = targ(2,:).' - src(2,:);
    r2 = rx.^2 + ry.^2;
    
    dx1 = repmat((targtang(1,:)).',1,ns);
    dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    tauxtarg = dx1./ds1;
    tauytarg = dy1./ds1;
    
    dx1 = repmat(srctang(1,:),nt,1);
    dy1 = repmat(srctang(2,:),nt,1);
    
    d2x1 = repmat(srcd2(1,:),nt,1);
    d2y1 = repmat(srcd2(2,:),nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    kappa = numer./denom;   
    
    xs = repmat(src(1,:),nt,1);
    ys = repmat(src(2,:),nt,1);
    
    % Starfish :  

    narms = 3;
    amp = 0.3;
    phi = pi/3;

    t = atan2(ys,xs);

    ct = cos(t);
    st = sin(t);
    cnt = cos(narms*(t + phi));
    snt = sin(narms*(t + phi));

    x1 = -(1+amp*cnt).*st-narms*amp*snt.*ct;
    y1 = (1+amp*cnt).*ct-narms*amp*snt.*st;

    x2 = -y1-narms*amp*(narms*cnt.*ct-snt.*st);
    y2 = x1-narms*amp*(narms*cnt.*st+snt.*ct);

    x3 = -y2 + amp*narms^3*snt.*ct + 2*amp*narms^2*cnt.*st + amp*narms*snt.*ct;
    y3 = x2 + amp*narms^3*snt.*st - 2*amp*narms^2*cnt.*ct + amp*narms*snt.*st;

    N = x1.*y2 - y1.*x2;
    dNdt = x1.*y3 - y1.*x3;
    D = (x1.^2 + y1.^2).^( 3/2);
    dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
    kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2); % kappa prime
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);      
    
    K11 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);
    
    K12 = -1/(2*zk^2).*(grad(:, :, 1).*nx + grad(:, :, 2).*ny);
           
    K22 = -1/(2*zk^2).*(third(:,:,1).*(nx.*nxtarg.^2) + third(:,:,2).*(2*nx.*nxtarg.*nytarg + ny.*nxtarg.^2) + ...
         + third(:,:,3).*(nx.*nytarg.^2 + 2*ny.*nxtarg.*nytarg) + third(:,:,4).*(ny.*nytarg.^2))+...
         -nu./(2*zk^2).*(third(:,:,1).*(nx.*tauxtarg.^2) + third(:,:,2).*(2*nx.*tauxtarg.*tauytarg + ny.*tauxtarg.^2) + ...
         + third(:,:,3).*(nx.*tauytarg.^2 + 2*ny.*tauxtarg.*tauytarg) + third(:,:,4).*(ny.*tauytarg.^2)); 
    
    
    [~, ~, ~, third, forth, fifth] = flex2d.hkdiffgreen(zk, src, targ, false);      
    
    K21 = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ... % K21 with the biharmonic Green's function subtracted off (biharmonic part handled in a separate function)
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}
    
    submat = zeros(2*nt,2*ns);
    
    submat(1:2:end,1:2:end) = K11;
    submat(1:2:end,2:2:end) = K12;
    
    submat(2:2:end,1:2:end) = K21; 
    submat(2:2:end,2:2:end) = K22; 

end

% supported plate kernel for modified biharmonic problem on a disk
if strcmpi(type, 'supported plate circle')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   srcd2 = srcinfo.d2;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};
   nu = coefs(1);

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   zkimag = (1i)*zk;

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   dx1 = repmat(srctang(1,:),nt,1);
   dy1 = repmat(srctang(2,:),nt,1);
   
   d2x1 = repmat(srcd2(1,:),nt,1);
   d2y1 = repmat(srcd2(2,:),nt,1);

   kappa = 1;   

   kp = 0;

   a1 = 2-nu;
   a2 = (-1+nu)*(7+nu)/(3 - nu);
   a3 = (1-nu)*(3+nu)/(1+nu);
   
   [~, grad, hess, third, ~] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part

   K11 = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

   K21 = 0;

   K12 = -1/(2*zk^2).*(grad(:, :, 1).*nx + grad(:, :, 2).*ny);
           
   K22 = -1/(2*zk^2).*(third(:,:,1).*(nx.*nxtarg.^2) + third(:,:,2).*(2*nx.*nxtarg.*nytarg + ny.*nxtarg.^2) + ...
         + third(:,:,3).*(nx.*nytarg.^2 + 2*ny.*nxtarg.*nytarg) + third(:,:,4).*(ny.*nytarg.^2))+...
         -nu./(2*zk^2).*(third(:,:,1).*(nx.*tauxtarg.^2) + third(:,:,2).*(2*nx.*tauxtarg.*tauytarg + ny.*tauxtarg.^2) + ...
         + third(:,:,3).*(nx.*tauytarg.^2 + 2*ny.*tauxtarg.*tauytarg) + third(:,:,4).*(ny.*tauytarg.^2)); 

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21; 
  submat(2:2:end,2:2:end) = K22; 

end

% supported plate kernel K21 for the flexural wave equation
if strcmpi(type, 'supported plate K21')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   srcd2 = srcinfo.d2;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};
   nu = coefs(1);

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   dx1 = repmat(srctang(1,:),nt,1);
   dy1 = repmat(srctang(2,:),nt,1);
   
   d2x1 = repmat(srcd2(1,:),nt,1);
   d2y1 = repmat(srcd2(2,:),nt,1);
   
   denom = sqrt(dx1.^2+dy1.^2).^3;
   numer = dx1.*d2y1-d2x1.*dy1;

   kappa = numer./denom;   

   xs = repmat(src(1,:),nt,1);
   ys = repmat(src(2,:),nt,1);

   % Ellipse:  

   c = 1;

   t = atan2(c*ys,xs); % determining 'theta' from the locations of x and y

   xs = c*cos(t);
   ys = sin(t);

   x1 = -c*sin(t);    
   y1 = cos(t);

   x2 = -c*cos(t);
   y2 = -sin(t);

   x3 = c*sin(t);
   y3 = -cos(t);

   % % Starfish :  
   % 
   % narms = 3;
   % amp = 0.3;
   % phi = pi/3;
   % 
   % t = atan2(ys,xs);
   % 
   % ct = cos(t);
   % st = sin(t);
   % cnt = cos(narms*(t + phi));
   % snt = sin(narms*(t + phi));
   % 
   % x1 = -(1+amp*cnt).*st-narms*amp*snt.*ct;
   % y1 = (1+amp*cnt).*ct-narms*amp*snt.*st;
   % 
   % x2 = -y1-narms*amp*(narms*cnt.*ct-snt.*st);
   % y2 = x1-narms*amp*(narms*cnt.*st+snt.*ct);
   % 
   % x3 = -y2 + amp*narms^3*snt.*ct + 2*amp*narms^2*cnt.*st + amp*narms*snt.*ct;
   % y3 = x2 + amp*narms^3*snt.*st - 2*amp*narms^2*cnt.*ct + amp*narms*snt.*st;
    
   N = x1.*y2 - y1.*x2;
   dNdt = x1.*y3 - y1.*x3;
   D = (x1.^2 + y1.^2).^( 3/2);
   dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
   kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2); % kappa prime

   a1 = 2-nu;
   a2 = (-1+nu)*(7+nu)/(3 - nu);
   a3 = (1-nu)*(3+nu)/(1+nu);
   
   [~, ~, ~, third, forth, ~] = flex2d.hkdiffgreen(zk, src, targ, false);      

   [~, ~, ~, ~, ~, fifth] = flex2d.hkdiffgreen(zk, src, targ, false);           

   submat = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}

end


% supported plate kernel K21 for the flexural wave equation
if strcmpi(type, 'supported plate K21 ellipse')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   srcd2 = srcinfo.d2;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};
   nu = coefs(1);

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   dx1 = repmat(srctang(1,:),nt,1);
   dy1 = repmat(srctang(2,:),nt,1);
   
   d2x1 = repmat(srcd2(1,:),nt,1);
   d2y1 = repmat(srcd2(2,:),nt,1);
   
   denom = sqrt(dx1.^2+dy1.^2).^3;
   numer = dx1.*d2y1-d2x1.*dy1;

   kappa = numer./denom;   

   xs = repmat(src(1,:),nt,1);
   ys = repmat(src(2,:),nt,1);

   % Ellipse:  

   c = 2;

   t = atan2(c*ys,xs); % determining 'theta' from the locations of x and y

   xs = c*cos(t);
   ys = sin(t);

   x1 = -c*sin(t);    
   y1 = cos(t);

   x2 = -c*cos(t);
   y2 = -sin(t);

   x3 = c*sin(t);
   y3 = -cos(t);

   N = x1.*y2 - y1.*x2;
   dNdt = x1.*y3 - y1.*x3;
   D = (x1.^2 + y1.^2).^( 3/2);
   dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
   kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2); % kappa prime

   a1 = 2-nu;
   a2 = (-1+nu)*(7+nu)/(3 - nu);
   a3 = (1-nu)*(3+nu)/(1+nu);
   
   % [~, ~, ~, third, forth, fifth] = flex2d.hkdiffgreen(zk, src, targ, true);      

   [~, ~, ~, thirdbh, forthbh, fifthbh] = flex2d.bhgreen(src,targ);
   third = 2*zk^2*thirdbh;
   forth = 2*zk^2*forthbh;
   fifth = 2*zk^2*fifthbh;

   submat = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}

end


% supported plate kernel K21 for the flexural wave equation
if strcmpi(type, 'supported plate K21 droplet')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    targnorm = targinfo.n;
    targtang = targinfo.d;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    rx = targ(1,:).' - src(1,:);
    ry = targ(2,:).' - src(2,:);
    r2 = rx.^2 + ry.^2;
    
    dx1 = repmat((targtang(1,:)).',1,ns);
    dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    tauxtarg = dx1./ds1;
    tauytarg = dy1./ds1;
    
    dx1 = repmat(srctang(1,:),nt,1);
    dy1 = repmat(srctang(2,:),nt,1);
    
    d2x1 = repmat(srcd2(1,:),nt,1);
    d2y1 = repmat(srcd2(2,:),nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    kappa = numer./denom;   
    
    xs = repmat(src(1,:),nt,1);
    ys = repmat(src(2,:),nt,1);
    
    % droplet: 
    
      a = 2;
    b = 1;
    c = -0.4;
    
    t = atan2(a*ys - c*xs.^2/a, b*xs);
    
    % x = a*cos(t);
    % y = b*sin(t) + c*cos(t).^2;
    
    x1 = -a*sin(t);
    y1 = b*cos(t) - 2*c*cos(t).*sin(t);
    
    x2 = -a*cos(t);
    y2 = -b*sin(t) + 2*c*sin(t).^2 - 2*c*cos(t).^2;
    
    x3 = a*sin(t);
    y3 = -b*cos(t) + 8*c*sin(t).*cos(t);
    
    N = x1.*y2 - y1.*x2;
    dNdt = x1.*y3 - y1.*x3;
    D = (x1.^2 + y1.^2).^( 3/2);
    dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
    kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2); % kappa prime
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    % [~, ~, ~, third, forth, fifth] = flex2d.hkdiffgreen(zk, src, targ, true);      
    
    [~, ~, ~, thirdbh, forthbh, fifthbh] = flex2d.bhgreen(src,targ);
    third = 2*zk^2*thirdbh;
    forth = 2*zk^2*forthbh;
    fifth = 2*zk^2*fifthbh;
    
    submat = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}

end



% supported plate kernel K21 for the flexural wave equation on a disk
if strcmpi(type, 'supported plate K21 circle')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   srcd2 = srcinfo.d2;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   coefs = varargin{1};
   nu = coefs(1);

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   dx1 = repmat(srctang(1,:),nt,1);
   dy1 = repmat(srctang(2,:),nt,1);
   
   d2x1 = repmat(srcd2(1,:),nt,1);
   d2y1 = repmat(srcd2(2,:),nt,1);

   kappa = 1;   
   kp = 0;

   a1 = 2-nu;
   a2 = (-1+nu)*(7+nu)/(3 - nu);
   a3 = (1-nu)*(3+nu)/(1+nu);
   
   [~, ~, ~, third, forth, ~] = flex2d.hkdiffgreen(zk, src, targ, false);      

   [~, ~, ~, ~, ~, fifth] = flex2d.hkdiffgreen(zk, src, targ, false);           

   submat = -1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*nxtarg.*nytarg.*nx.^3+3*nxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(nytarg.^2.*nx.^3 + 6*nxtarg.*nytarg.*nx.^2.*ny+3.*nxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*nytarg.^2.*nx.^2.*ny+6.*nxtarg.*nytarg.*nx.*ny.^2+nxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(nytarg.*ny.^2.*(3*nytarg.*nx+2*nxtarg.*ny))+ ...
          fifth(:,:,6).*nytarg.^2.*ny.^3) + ... % G_{nx nx ny ny ny}
         -nu./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.^3 + ...
          fifth(:,:,2).*(2*tauxtarg.*tauytarg.*nx.^3+3*tauxtarg.^2.*nx.^2.*ny) + ...
          fifth(:,:,3).*(tauytarg.^2.*nx.^3 + 6*tauxtarg.*tauytarg.*nx.^2.*ny+3.*tauxtarg.^2.*nx.*ny.^2) + ...
          fifth(:,:,4).*(3*tauytarg.^2.*nx.^2.*ny+6.*tauxtarg.*tauytarg.*nx.*ny.^2+tauxtarg.^2.*ny.^3)+ ...
          fifth(:,:,5).*(tauytarg.*ny.^2.*(3*tauytarg.*nx+2*tauxtarg.*ny))+ ...
          fifth(:,:,6).*tauytarg.^2.*ny.^3) + ... % nu*G_{taux taux ny ny ny}
        -a1./(2*zk^2).*(fifth(:,:,1).*nxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(2.*nx.*nxtarg.^2.*taux.*tauy + 2.*nx.*nxtarg.*nytarg.*taux.^2 + nxtarg.^2.*ny.*taux.^2) + ...
          fifth(:,:,3).*(nx.*nxtarg.^2.*tauy.^2 + 4*nx.*nxtarg.*nytarg.*taux.*tauy + 2*nxtarg.^2.*ny.*taux.*tauy + nx.*nytarg.^2.*taux.^2 + 2*nxtarg.*ny.*nytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*nxtarg.*nytarg.*tauy.^2 + nxtarg.^2.*ny.*tauy.^2 + 2*nx.*nytarg.^2.*taux.*tauy + 4*nxtarg.*ny.*nytarg.*taux.*tauy + ny.*nytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(nytarg.*tauy.*(2*ny.*nytarg.*taux + 2*nxtarg.*ny.*tauy + nx.*nytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*nytarg.^2.*tauy.^2) + ...  % G_{nx nx ny tauy tauy}
        -nu*a1./(2*zk^2).*(fifth(:,:,1).*tauxtarg.^2.*nx.*taux.^2 + ...
          fifth(:,:,2).*(ny.*taux.^2.*tauxtarg.^2 + 2.*nx.*taux.*tauxtarg.^2.*tauy +2.*nx.*taux.^2.*tauxtarg.*tauytarg) + ...
          fifth(:,:,3).*(nx.*tauxtarg.^2.*tauy.^2 + 4*nx.*tauxtarg.*tauytarg.*taux.*tauy + 2*tauxtarg.^2.*ny.*taux.*tauy + nx.*tauytarg.^2.*taux.^2 + 2*tauxtarg.*ny.*tauytarg.*taux.^2) + ...
          fifth(:,:,4).*(2*nx.*tauxtarg.*tauytarg.*tauy.^2 + tauxtarg.^2.*ny.*tauy.^2 + 2*nx.*tauytarg.^2.*taux.*tauy + 4*tauxtarg.*ny.*tauytarg.*taux.*tauy + ny.*tauytarg.^2.*taux.^2) + ...
          fifth(:,:,5).*(tauytarg.*tauy.*(2*ny.*tauytarg.*taux + 2*tauxtarg.*ny.*tauy + nx.*tauytarg.*tauy)) + ...
          fifth(:,:,6).*ny.*tauytarg.^2.*tauy.^2) + ... % nu*G_{taux taux ny tauy tauy}
         a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(nxtarg.*nxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*nxtarg.*nytarg + 2*nx.*nxtarg.^2.*ny) + ...
          forth(:, :, 3).*(nxtarg.^2.*ny.^2 + 4*nx.*nxtarg.*ny.*nytarg + nx.^2.*nytarg.^2) +...
          forth(:, :, 4).*(2*ny.*nytarg.*(nxtarg.*ny + nx.*nytarg))+...
          forth(:, :, 5).*(nytarg.*nytarg.*ny.*ny)) + ... % kappa*G_{nx nx ny ny}
         nu*a2.*kappa./(2*zk^2).*(forth(:, :, 1).*(tauxtarg.*tauxtarg.*nx.*nx) + ...
          forth(:, :, 2).*(2*nx.^2.*tauxtarg.*tauytarg + 2*nx.*tauxtarg.^2.*ny) + ...
          forth(:, :, 3).*(tauxtarg.^2.*ny.^2 + 4*nx.*tauxtarg.*ny.*tauytarg + nx.^2.*tauytarg.^2) +...
          forth(:, :, 4).*(2*ny.*tauytarg.*(tauxtarg.*ny + nx.*tauytarg))+...
          forth(:, :, 5).*(tauytarg.*tauytarg.*ny.*ny)) + ... % kappa*nu*G_{nx nx ny ny}
        -a3.*kp./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
            third(:, :, 4).*(nytarg.*nytarg.*tauy)) + ... % kp*G_{nx nx tauy}
         -nu*a3.*kp./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*tauy)) ;  % kp*nu*G_{taux taux tauy}

end


% first kernel of evaluation for the support plate
if strcmpi(type, 'supported plate K1 eval')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   srcd2 = srcinfo.d2;
   coefs = varargin{1};
   nu = coefs(1);

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   

   zkimag = (1i)*zk;

   [~,grad] = chnk.lap2d.green(src,targ,true);
   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx); 

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   dx1 = repmat(srctang(1,:),nt,1);
   dy1 = repmat(srctang(2,:),nt,1);
   
   d2x1 = repmat(srcd2(1,:),nt,1);
   d2y1 = repmat(srcd2(2,:),nt,1);
   
   denom = sqrt(dx1.^2+dy1.^2).^3;
   numer = dx1.*d2y1-d2x1.*dy1;

   kappa = numer./denom; 

   xs = repmat(src(1,:),nt,1);
   ys = repmat(src(2,:),nt,1);

   % Ellipse:
   c = 1;

   t = atan2(c*ys,xs);

   x1 = -c*sin(t);    
   y1 = cos(t);

   x2 = -c*cos(t);
   y2 = -sin(t);

   x3 = c*sin(t);
   y3 = -cos(t);

   % % Starfish :  
   % 
   % narms = 3;
   % amp = 0.3;
   % phi = pi/3;
   % 
   % t = atan2(ys,xs);
   % 
   % ct = cos(t);
   % st = sin(t);
   % cnt = cos(narms*(t + phi));
   % snt = sin(narms*(t + phi));
   % 
   % x1 = -(1+amp*cnt).*st-narms*amp*snt.*ct;
   % y1 = (1+amp*cnt).*ct-narms*amp*snt.*st;
   % 
   % x2 = -y1-narms*amp*(narms*cnt.*ct-snt.*st);
   % y2 = x1-narms*amp*(narms*cnt.*st+snt.*ct);
   % 
   % x3 = -y2 + amp*narms^3*snt.*ct + 2*amp*narms^2*cnt.*st + amp*narms*snt.*ct;
   % y3 = x2 + amp*narms^3*snt.*st - 2*amp*narms^2*cnt.*ct + amp*narms*snt.*st;

   N = x1.*y2 - y1.*x2;
   dNdt = x1.*y3 - y1.*x3;
   D = (x1.^2 + y1.^2).^( 3/2);
   dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
   kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2);

   a1 = 2-nu;
   a2 = (-1+nu)*(7+nu)/(3 - nu);
   a3 = (1-nu)*(3+nu)/(1+nu);

   [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part


   %[~, grad, hess, third, ~, ~] = flex2d.helmdiffgreen(zk, src, targ, false);            % Hankel part
   %[~, gradK, hessK, thirdK, ~, ~] = flex2d.helmdiffgreen(zkimag, src, targ, false);     % modified bessel K part 
   

   submat = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

end



% first kernel of evaluation for the support plate
if strcmpi(type, 'supported plate K1 eval')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   srcd2 = srcinfo.d2;
   coefs = varargin{1};
   nu = coefs(1);

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   

   zkimag = (1i)*zk;

   [~,grad] = chnk.lap2d.green(src,targ,true);
   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx); 

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   dx1 = repmat(srctang(1,:),nt,1);
   dy1 = repmat(srctang(2,:),nt,1);
   
   d2x1 = repmat(srcd2(1,:),nt,1);
   d2y1 = repmat(srcd2(2,:),nt,1);
   
   denom = sqrt(dx1.^2+dy1.^2).^3;
   numer = dx1.*d2y1-d2x1.*dy1;

   kappa = numer./denom; 

   xs = repmat(src(1,:),nt,1);
   ys = repmat(src(2,:),nt,1);

   % Ellipse:
   c = 1;

   t = atan2(c*ys,xs);

   x1 = -c*sin(t);    
   y1 = cos(t);

   x2 = -c*cos(t);
   y2 = -sin(t);

   x3 = c*sin(t);
   y3 = -cos(t);

   % % Starfish :  
   % 
   % narms = 3;
   % amp = 0.3;
   % phi = pi/3;
   % 
   % t = atan2(ys,xs);
   % 
   % ct = cos(t);
   % st = sin(t);
   % cnt = cos(narms*(t + phi));
   % snt = sin(narms*(t + phi));
   % 
   % x1 = -(1+amp*cnt).*st-narms*amp*snt.*ct;
   % y1 = (1+amp*cnt).*ct-narms*amp*snt.*st;
   % 
   % x2 = -y1-narms*amp*(narms*cnt.*ct-snt.*st);
   % y2 = x1-narms*amp*(narms*cnt.*st+snt.*ct);
   % 
   % x3 = -y2 + amp*narms^3*snt.*ct + 2*amp*narms^2*cnt.*st + amp*narms*snt.*ct;
   % y3 = x2 + amp*narms^3*snt.*st - 2*amp*narms^2*cnt.*ct + amp*narms*snt.*st;

   N = x1.*y2 - y1.*x2;
   dNdt = x1.*y3 - y1.*x3;
   D = (x1.^2 + y1.^2).^( 3/2);
   dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
   kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2);

   a1 = 2-nu;
   a2 = (-1+nu)*(7+nu)/(3 - nu);
   a3 = (1-nu)*(3+nu)/(1+nu);

   [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part


   %[~, grad, hess, third, ~, ~] = flex2d.helmdiffgreen(zk, src, targ, false);            % Hankel part
   %[~, gradK, hessK, thirdK, ~, ~] = flex2d.helmdiffgreen(zkimag, src, targ, false);     % modified bessel K part 
   

   submat = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

end

% first kernel of evaluation for the support plate
if strcmpi(type, 'supported plate K1 eval ellipse')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   srcd2 = srcinfo.d2;
   coefs = varargin{1};
   nu = coefs(1);

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   

   zkimag = (1i)*zk;

   [~,grad] = chnk.lap2d.green(src,targ,true);
   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx); 

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   dx1 = repmat(srctang(1,:),nt,1);
   dy1 = repmat(srctang(2,:),nt,1);
   
   d2x1 = repmat(srcd2(1,:),nt,1);
   d2y1 = repmat(srcd2(2,:),nt,1);
   
   denom = sqrt(dx1.^2+dy1.^2).^3;
   numer = dx1.*d2y1-d2x1.*dy1;

   kappa = numer./denom; 

   xs = repmat(src(1,:),nt,1);
   ys = repmat(src(2,:),nt,1);

   % Ellipse:
   c = 2;

   t = atan2(c*ys,xs);

   x1 = -c*sin(t);    
   y1 = cos(t);

   x2 = -c*cos(t);
   y2 = -sin(t);

   x3 = c*sin(t);
   y3 = -cos(t);

   N = x1.*y2 - y1.*x2;
   dNdt = x1.*y3 - y1.*x3;
   D = (x1.^2 + y1.^2).^( 3/2);
   dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
   kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2);

   a1 = 2-nu;
   a2 = (-1+nu)*(7+nu)/(3 - nu);
   a3 = (1-nu)*(3+nu)/(1+nu);

   [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part


   submat = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

end


% first kernel of evaluation for the support plate
if strcmpi(type, 'supported plate K1 eval droplet')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    
    
    zkimag = (1i)*zk;
    
    [~,grad] = chnk.lap2d.green(src,targ,true);
    hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx); 
    
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    rx = targ(1,:).' - src(1,:);
    ry = targ(2,:).' - src(2,:);
    r2 = rx.^2 + ry.^2;
    
    dx1 = repmat(srctang(1,:),nt,1);
    dy1 = repmat(srctang(2,:),nt,1);
    
    d2x1 = repmat(srcd2(1,:),nt,1);
    d2y1 = repmat(srcd2(2,:),nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    kappa = numer./denom; 
    
    xs = repmat(src(1,:),nt,1);
    ys = repmat(src(2,:),nt,1);
    
    % droplet: 
    
    a = 2;
    b = 1;
    c = -0.4;
    
    t = atan2(a*ys - c*xs.^2/a, b*xs);
    
    % x = a*cos(t);
    % y = b*sin(t) + c*cos(t).^2;
    
    x1 = -a*sin(t);
    y1 = b*cos(t) - 2*c*cos(t).*sin(t);
    
    x2 = -a*cos(t);
    y2 = -b*sin(t) + 2*c*sin(t).^2 - 2*c*cos(t).^2;
    
    x3 = a*sin(t);
    y3 = -b*cos(t) + 8*c*sin(t).*cos(t);
    
    N = x1.*y2 - y1.*x2;
    dNdt = x1.*y3 - y1.*x3;
    D = (x1.^2 + y1.^2).^( 3/2);
    dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
    kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2);
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part
    
    
    submat = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

end

% first kernel of evaluation for the support plate
if strcmpi(type, 'supported plate K1 eval starfish 3arms 1')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    
    
    zkimag = (1i)*zk;
    
    [~,grad] = chnk.lap2d.green(src,targ,true);
    hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx); 
    
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    rx = targ(1,:).' - src(1,:);
    ry = targ(2,:).' - src(2,:);
    r2 = rx.^2 + ry.^2;
    
    dx1 = repmat(srctang(1,:),nt,1);
    dy1 = repmat(srctang(2,:),nt,1);
    
    d2x1 = repmat(srcd2(1,:),nt,1);
    d2y1 = repmat(srcd2(2,:),nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    kappa = numer./denom; 
    
    xs = repmat(src(1,:),nt,1);
    ys = repmat(src(2,:),nt,1);
    
    % Starfish :  

    narms = 3;
    amp = 0.3;
    phi = 0;

    t = atan2(ys,xs);

    ct = cos(t);
    st = sin(t);
    cnt = cos(narms*(t + phi));
    snt = sin(narms*(t + phi));

    x1 = -(1+amp*cnt).*st-narms*amp*snt.*ct;
    y1 = (1+amp*cnt).*ct-narms*amp*snt.*st;

    x2 = -y1-narms*amp*(narms*cnt.*ct-snt.*st);
    y2 = x1-narms*amp*(narms*cnt.*st+snt.*ct);

    x3 = -y2 + amp*narms^3*snt.*ct + 2*amp*narms^2*cnt.*st + amp*narms*snt.*ct;
    y3 = x2 + amp*narms^3*snt.*st - 2*amp*narms^2*cnt.*ct + amp*narms*snt.*st;
    
    N = x1.*y2 - y1.*x2;
    dNdt = x1.*y3 - y1.*x3;
    D = (x1.^2 + y1.^2).^( 3/2);
    dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
    kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2);
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part
    
    
    submat = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

end


% first kernel of evaluation for the support plate
if strcmpi(type, 'supported plate K1 eval starfish 3arms 2')
    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    srcd2 = srcinfo.d2;
    coefs = varargin{1};
    nu = coefs(1);
    
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    
    
    zkimag = (1i)*zk;
    
    [~,grad] = chnk.lap2d.green(src,targ,true);
    hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx); 
    
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    
    ds = sqrt(dx.*dx+dy.*dy);
    
    taux = dx./ds;
    tauy = dy./ds;
    
    rx = targ(1,:).' - src(1,:);
    ry = targ(2,:).' - src(2,:);
    r2 = rx.^2 + ry.^2;
    
    dx1 = repmat(srctang(1,:),nt,1);
    dy1 = repmat(srctang(2,:),nt,1);
    
    d2x1 = repmat(srcd2(1,:),nt,1);
    d2y1 = repmat(srcd2(2,:),nt,1);
    
    denom = sqrt(dx1.^2+dy1.^2).^3;
    numer = dx1.*d2y1-d2x1.*dy1;
    
    kappa = numer./denom; 
    
    xs = repmat(src(1,:),nt,1);
    ys = repmat(src(2,:),nt,1);
    
    % Starfish :  

    narms = 3;
    amp = 0.3;
    phi = pi/3;

    t = atan2(ys,xs);

    ct = cos(t);
    st = sin(t);
    cnt = cos(narms*(t + phi));
    snt = sin(narms*(t + phi));

    x1 = -(1+amp*cnt).*st-narms*amp*snt.*ct;
    y1 = (1+amp*cnt).*ct-narms*amp*snt.*st;

    x2 = -y1-narms*amp*(narms*cnt.*ct-snt.*st);
    y2 = x1-narms*amp*(narms*cnt.*st+snt.*ct);

    x3 = -y2 + amp*narms^3*snt.*ct + 2*amp*narms^2*cnt.*st + amp*narms*snt.*ct;
    y3 = x2 + amp*narms^3*snt.*st - 2*amp*narms^2*cnt.*ct + amp*narms*snt.*st;
    
    N = x1.*y2 - y1.*x2;
    dNdt = x1.*y3 - y1.*x3;
    D = (x1.^2 + y1.^2).^( 3/2);
    dDdt = 3/2.*((x1.^2 + y1.^2).^(1/2)).*(2.*x1.*x2 + 2.*y1.*y2);
    kp = ( D.*dNdt - N.*dDdt) ./ (x1.^2 + y1.^2).^( 7/2);
    
    a1 = 2-nu;
    a2 = (-1+nu)*(7+nu)/(3 - nu);
    a3 = (1-nu)*(3+nu)/(1+nu);
    
    [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part
    
    
    submat = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

end


% first kernel of evaluation for the support plate
if strcmpi(type, 'supported plate K1 eval circle')
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   srcd2 = srcinfo.d2;
   coefs = varargin{1};
   nu = coefs(1);

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   

   zkimag = (1i)*zk;

   [~,grad] = chnk.lap2d.green(src,targ,true);
   hilb = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx); 

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;

   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   dx1 = repmat(srctang(1,:),nt,1);
   dy1 = repmat(srctang(2,:),nt,1);
   
   d2x1 = repmat(srcd2(1,:),nt,1);
   d2y1 = repmat(srcd2(2,:),nt,1);

   kappa = 1;
   kp = 0;

   a1 = 2-nu;
   a2 = (-1+nu)*(7+nu)/(3 - nu);
   a3 = (1-nu)*(3+nu)/(1+nu);

   [~, grad, hess, third, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part

   submat = -1/(2*zk^2)*(third(:,:,1).*nx.^3 + 3*third(:,:,2).*nx.^2.*ny + 3*third(:,:,3).*nx.*ny.^2 + third(:,:,4).*ny.^3) + ...
         -a1/(2*zk^2)*(third(:,:,1).*nx.*taux.^2 + third(:,:,2).*(ny.*taux.^2 + 2*nx.*taux.*tauy) + third(:,:,3).*(nx.*tauy.^2 + 2*ny.*taux.*tauy) + third(:,:,4).*ny.*tauy.^2) + ...
         a2*kappa./(2*zk^2).*(hess(:,:,1).*nx.^2 + 2*hess(:,:,2).*nx.*ny + hess(:,:,3).*ny.^2) + ...
         -a3*kp./(2*zk^2).*(grad(:,:,1).*taux + grad(:,:,2).*tauy);

end



% second kernel of evaluation for the supported plate
if strcmpi(type, 'supported plate K2 eval')
    srcnorm = srcinfo.n;

    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);

    
    [~, grad,~, ~, ~, ~] = flex2d.hkdiffgreen(zk, src, targ, false);            % Hankel part
    

    submat = 1/(2*zk^2).*(grad(:,:,1).*nx + grad(:,:,2).*ny);

    submat = -submat;
end


