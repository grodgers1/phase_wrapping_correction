%% name: dpcprocess
%%
%% syntax: [amp,dpc,dci] = dpcprocess(ffm2,datm2,nperiods,mf2);
%%
%% description: dpc image processing
%%
%%		ffm2,datm2 - input/output 3D data set (phasesteps along 3 dimension)
%%		nperiods - number of oscillations in period
%%		mf2 - median filter after processing ([1 1] or [3 3])

function [amp,dpc,dci,visff,visdat] = dpcprocess_griffin(ffm2,datm2,nperiods,mf2,dowrap);
if (nargin < 4)
    fprintf('Usage:\n');
    fprintf('[amp,dpc,dci,visff,visdat] = dpcprocess(ffm2,datm2,nperiods,mf2);\n');
    return;
end

%% Fourier transform & extract angles

fftffm2 = fft(ffm2(:,:,:),[],3); clear ffm2;
fftdatm2 = fft(datm2(:,:,:),[],3); clear datm2;

ampdat = abs(fftdatm2(:,:,1)); 
ampff = abs(fftffm2(:,:,1)); 
dpcdat = angle(fftdatm2(:,:,1+nperiods)); 
dpcff = angle(fftffm2(:,:,1+nperiods)); 
dcidat = abs(fftdatm2(:,:,1+nperiods));
dciff = abs(fftffm2(:,:,1+nperiods));
visff = abs(fftffm2(:,:,1+nperiods))./abs(fftffm2(:,:,1));
visdat = abs(fftdatm2(:,:,1+nperiods))./abs(fftdatm2(:,:,1));

clear fftffm2; clear fftdatm2;

%% filter again & wrap difference data

if mf2 == 1
    amp = ampdat./ampff; 
	if dowrap == 1
        dpc = wrap(dpcdat-dpcff,2);
    else
        dpc = dpcdat-dpcff; 
    end
    dci = (dcidat./ampdat)./(dciff./ampff);
else
    amp = medfilt2(ampdat./ampff,[mf2 mf2]); 
	if dowrap == 1
        dpc = medfilt2(wrap(dpcdat-dpcff,2),[mf2 mf2]); 
    else
        dpc = medfilt2(dpcdat-dpcff,[mf2 mf2]); 
    end
	dci = medfilt2((dcidat./ampdat)./(dciff./ampff),[mf2 mf2]);
end
