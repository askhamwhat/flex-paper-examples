function [sysmat] = get_trans_mat(src,tar,kerns)

    nt = size(tar.r,2);
    ns = size(src.r,2);

    fkern = kerns{1};
    fkern1 = kerns{2};
    fkern2 = kerns{3};
    hilbert = kerns{4};
    double = kerns{5};
    fkern3 = kerns{6};
    fkern4 = kerns{7};
    fkern5 = kerns{8};
    fkern6 = kerns{9};
    fkern1_trans = kerns{10};
    fkern3_trans = kerns{11};

    ikern1 = kerns{12};
    ikern2 = kerns{13};
    ikern3 = kerns{14};

    wts = src.wts(:).';

    sysmat = zeros(nt*2,ns*3);

    fsys = fkern(src,tar);
    sysmat(:,1:3:end) = fsys(:,1:2:end).*wts;
    sysmat(:,3:3:end) = fsys(:,2:2:end).*wts;

    sysmat1 = fkern1_trans(src,tar).*wts;
    sysmat2 = fkern2(src,tar).*wts;

    K21 = fkern3_trans(src,tar).*wts;
    K21second = fkern4(src,tar).*wts;
    K21hilbert = fkern5(src,tar).*wts;
    K22second = fkern6(src,tar).*wts;

    kappa = tar.kappa;
    %kappa = signed_curvature(tar);
    %kappa = kappa(:)';

    hilb = sysmat1;
    hilb2 = sysmat2;

    mat1 = sysmat(1:2:end, 1:3:end);
    mat4 = sysmat(2:2:end, 3:3:end);

    sysmat(1:2:end, 1:3:end) = mat1;
    sysmat(2:2:end, 1:3:end) = K21  + (kappa.').*(K21second);
    sysmat(2:2:end, 3:3:end) = mat4 + (kappa.').*(K22second);

    sysmat(1:2:end,2:3:end) = hilb;
    sysmat(2:2:end,2:3:end) = hilb2 + (kappa.').*(K21hilbert);

end