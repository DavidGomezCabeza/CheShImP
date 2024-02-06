function entr = phaseCorrectCostFunction(spectrum, params)

  Rpre = real(spectrum);

  PeakInfo  = GetPeaks( Rpre, 6, 0.001 );
  Weight = ones( length( spectrum),1 );
  for i = 1 : length( PeakInfo )
    Weight( PeakInfo( i ).Start:PeakInfo( i ).End ) = 0;
  end


  L = length(Rpre);
  E = speye( L );
  D = diff( E, 1 );
  W = spdiags(Weight, 0, L, L );
  C = chol( W + 1600 * D' * D );
  try
    EntropyBaseLine = C\( C'\( Weight.* Rpre' ) );
  catch
    EntropyBaseLine = 0;
  end
 % baseline correction
 R = Rpre - EntropyBaseLine';

  % % % calculate Normalized histogram (entropy)
  h = abs(R) / sum(abs(R)); % normalized histogram of the absolute values of the real part (probability distribution)
  entr  = -sum(h .* log(h)); % calculate the entropy 


  % apply non-negativity penalty
  
  if(params.nonNegativePenalty)
    % equation 8 in reference
    stepDownFunction = (0.5 - sign(R)/2); % 1 when negative, 0 when positive 

    % sign returns an array y the same size as x, where each element of y
    % is 1 if the corresponding element of x is greater than 0. 0 if the
    % corresponding element of x equals 0. and -1 if the corresponding
    % element of x is less than 0. if R is negative, we have 1 because sign
    % give out -1, so becames +1 with the minus before. 
  
    % should be determined to make negativity penalty on order of entropy
    % for now hard coded
    scaleFactor = 1e-1;

    % equation 7 in reference
    % should be abs(R) instead of abs(h)? better check normalization
    negativePenalty = scaleFactor * sum(stepDownFunction .* abs(R));
    entr = entr + negativePenalty;
  end
  % 
  if max(R) ~= max(abs(R))
    entr = entr*2;
  end

end
