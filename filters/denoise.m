function stftdn = denoise(stft,thresh)

ang = angle(stft); % Save phase

stftdb = mag2db(abs(stft)); % Convert signal amplitude module to dB

stftdbdn = (stftdb<thresh)*(-90) + stftdb; % Reduce dB under thershold by -90 dB
%sum(stftdb<thresh)
stftdn = db2mag(stftdbdn); % Get back to amplitude

stftdn = stftdn.*(cos(ang)+1i.*sin(ang)); % Get back to complex signal (with phase information)

end
