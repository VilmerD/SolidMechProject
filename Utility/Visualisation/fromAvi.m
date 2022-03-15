function fromAvi(name, format)
files = dir(sprintf('**/%s', name));
input_file_folder = files(1).folder;
input_file_name = files(1).name;
input_fullfile = fullfile(input_file_folder, input_file_name);

output_file_name = input_file_name(1:strfind(input_file_name, '.avi')-1);
output_fullfile = fullfile(input_file_folder, output_file_name);

input_reader = VideoReader(input_fullfile);
output_writer = VideoWriter(output_fullfile, format);
output_writer.FrameRate = input_reader.FrameRate;
open(output_writer);

for count = 1:input_reader.NumFrames
    frame = read(input_reader, count);
    writeVideo(output_writer, frame)
end

close(input_reader);
close(output_writer);
end