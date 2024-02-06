function filterPeak = GetPeaks(Spec, NoiseFactor, CwtFilterFactor)
%  Input 
%         AverageData: the real data of speactrum data
%         NoiseFactor: noise factor
%         CwtFilterFactor: cwt factor
%  Output
%         filterPeak: peak information
% 

% DA ELIMINARE -------
% Spec= DataBeforePhase;  
% nNoiseFactor = 6; 
% CwtFilterFactor= 0.0008; 
% --------------------

chooseFID = 0; 
StartOfNoise = 1; 

SpecData = real (Spec); % real part of the spectrum 
dataSize = length(SpecData); % dimension of the data 

% Continues Wavelet Transform 
RealDiffSpec = cwt( real(Spec), CwtFilterFactor * dataSize, 'haar' ); % diff because it's the first derivate 
ImagDiffSpec = cwt( imag(Spec), CwtFilterFactor * dataSize, 'haar' );

DiffData = abs(RealDiffSpec+sqrt(-1)*ImagDiffSpec); % sqrt(-1)=i
dataDiffSize = length(DiffData); 

% calcualte noise level 
if (chooseFID==0)
    NoisePoint = round(dataSize/16); 
    HalfNoisePoint=round(NoisePoint/2); 
    tempNoise = sum(DiffData(StartOfNoise:StartOfNoise+NoisePoint-1).^2);
    temp=(sum(DiffData(StartOfNoise:StartOfNoise+NoisePoint-1)))^2;
    tempSym=0; 
    for i = 1:HalfNoisePoint
        tempSym = tempSym + i  *  (DiffData(StartOfNoise + HalfNoisePoint + i-1) - DiffData(StartOfNoise + HalfNoisePoint - i));
    end
    tempSym = tempSym  *  3 / ( NoisePoint  *  NoisePoint - 1 );
    temp2 = temp + tempSym;
    temp2 = temp2 / NoisePoint;
    Noise = sqrt(( tempNoise - temp2 ) / ( NoisePoint - 1 ));
else % the second method to calculate the noise level divides the spectrogram in 16 parts that are the same and found the mean square value for each segment, taking the smaller value
    NoisePoint = round(dataDiffSize / 16);
    tempNoise = zeros( 1,16 );
    tempAvg = zeros( 1,16 );
    for i = 1:16
        tempNoise(1,i) = std( DiffData( StartOfNoise + (i-1) * NoisePoint : StartOfNoise +  (i)*NoisePoint - 1 ) );
        tempAvg(1,i) = mean( DiffData( StartOfNoise + (i - 1) * NoisePoint : StartOfNoise + (i)*NoisePoint - 1));
    end
    Noise = min( tempNoise );
end

% using sliding window to recognize the peak in spectra
Noise=Noise * NoiseFactor;
PeakPickState = 0; % 0 = means no peaks (baseline); 1 = means peak start point; 2 = menas peak end point reached
PeakMaxSlope = 0;
TempPeakStart = 0;
TempPeakEnd = 0;
TempPeakSlope = 0;
TempPeakIndex = 0;
PeakNum = 1; % maybe it should be changed
HalfWindow = round( CwtFilterFactor * dataSize/2 );

% below, for the sake of simplicity, the first 10 points (our case 7) are considered to
% be necessarily noisy, so we conder them directly from the 11th point  (in
% out case 8th point). The same for the last point

for i = ( HalfWindow +  1 ) : length( DiffData ) - HalfWindow
   MaxOfWindow = max( DiffData( i - HalfWindow : i +  HalfWindow ) ); % maximum value of the window 
   MinOfWindow = min( DiffData( i - HalfWindow : i +  HalfWindow ) ); % minimum value 
   WindowHeight = MaxOfWindow-MinOfWindow; 
   if PeakPickState == 0 % base line 
       if WindowHeight>Noise
           TempPeakStart = i + HalfWindow-5; 
           PeakPickState = 1; % peak is started
           if( PeakMaxSlope<DiffData( i ) )
                TempPeakIndex = i;
                TempPeakSlope = WindowHeight;
           end
       end
   end

   if ( ( PeakPickState == 1 ) )
        if( TempPeakSlope <  WindowHeight )
                TempPeakIndex = i;
                TempPeakSlope = WindowHeight;
        end
       if( i == dataDiffSize || WindowHeight < Noise )
           PeakPickState = 0;
           TempPeakEnd = i - HalfWindow +  5;
           PeakInfo( PeakNum ) = struct( 'Start', TempPeakStart, 'End', TempPeakEnd, 'Index', TempPeakIndex, 'MaxSlope', TempPeakSlope, 'Height', 0 );
           TempPeakSlope = 0;
           PeakNum = PeakNum +  1;
       end
   end
end
   % this is where the start and end of the peaks have been found 

filterPeak = [];
try
for i = 1 : length( PeakInfo )
    if( ( PeakInfo( i ).End - PeakInfo( i ).Start ) > 10 )
        filterPeak = [ filterPeak, PeakInfo( i ) ];
    end
end
catch
end

