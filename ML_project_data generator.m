n=50;
data=zeros(n,32768+4096);
for i=1:n
    
    % Magnetic induction
    sys.magnet=14.095;

    % Spin system
    sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'};

    % Chemical shifts
    inter.zeeman.scalar={10*rand(1),10*rand(1),10*rand(1),10*rand(1),10*rand(1),10*rand(1),10*rand(1),...
                        10*rand(1),10*rand(1),10*rand(1),10*rand(1),10*rand(1)};
    % 
    % Scalar couplings
    inter.coupling.scalar=cell(12,12);
    inter.coupling.scalar{1,2}=10*rand(1);
    inter.coupling.scalar{1,3}=10*rand(1);
    inter.coupling.scalar{2,3}=10*rand(1);
    inter.coupling.scalar{4,5}=10*rand(1);
    inter.coupling.scalar{4,6}=10*rand(1);
    inter.coupling.scalar{5,6}=10*rand(1);
    inter.coupling.scalar{8,9}=10*rand(1);
    inter.coupling.scalar{8,10}=10*rand(1);
    inter.coupling.scalar{9,10}=10*rand(1);
    inter.coupling.scalar{10,11}=10*rand(1);
    inter.coupling.scalar{10,12}=10*rand(1);
    inter.coupling.scalar{11,12}=10*rand(1);

    % Basis set
    bas.formalism='sphten-liouv';
    bas.approximation='IK-2';
    bas.connectivity='scalar_couplings';
    bas.space_level=1;

    % Spinach housekeeping
    spin_system=create(sys,inter);
    spin_system=basis(spin_system,bas);

    % Sequence parameters
    parameters.spins={'1H'};
    parameters.rho0=state(spin_system,'L+','1H');
    parameters.coil=state(spin_system,'L+','1H');
    parameters.decouple={};
    parameters.offset=4600;
    parameters.sweep=1200;
    parameters.npoints=4096;
    parameters.zerofill=32768;
    parameters.axis_units='ppm';
    parameters.invert_axis=1;

    % Simulation
    fid=liquid(spin_system,@acquire,parameters,'nmr');

    % Apodization
    fid=apodization(fid,'exp-1d',10);
    fid=awgn(fid,50);
    % Fourier transform
    spectrum=fftshift(fft(fid,parameters.zerofill));
    
    % sample of data
    data(i,1:4096)=fid;
    data(i,4097:36864)=spectrum;

end