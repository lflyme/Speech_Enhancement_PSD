function speech_est = overlap_add(signal)

[fr_size , fr_count] = size(signal);
fr_overlap = 0.5 * fr_size;

speech_est = zeros(length(1:fr_count*fr_overlap), 1);
counter = 1;
for k = 1 : fr_overlap : (fr_count - 1) * fr_overlap
    % Frame selection - Current index to next (fr_size) indices
    fr_unit = k : k + fr_size - 1;
    speech_est(fr_unit) = speech_est(fr_unit) + signal(:,counter);
    counter = counter + 1;
end
end