%% 帧差法
video = VideoReader('768x576.avi');
nFrames = video.NumberOfFrames;
for i = 1:nFrames-1
    frame_pre = read(video, i); %前一帧
    frame_pre = rgb2gray(frame_pre);
    frame_pos = read(video, i+1); %当前帧
    frame_pos = rgb2gray(frame_pos);
    current_object = abs(frame_pre - frame_pos); % 作差
    current_object = im2bw(current_object,0.2);
    imshow(current_object);
end

%% 单高斯背景建模
p = 0.9;
frame_temp = read(video, 1);
frame_temp = im2double(rgb2gray(frame_temp));
[m,n] = size(frame_temp);
% avg_mat = zeros(m, n, 3);
% var_mat = zeros(m, n, 3);
gray_100 = zeros(m,n,100);
for i = 1:100 % 取前100帧统计均值方差
    frame_temp = read(video, i);
    frame_temp = im2double(rgb2gray(frame_temp));
    gray_100(:,:,i) = frame_temp;
end

avg_mat = mean(gray_100, 3);
std_mat = std(gray_100, 0, 3);
clear gray_100
b_pre = avg_mat;

for i = 101:nFrames
    frame_temp = read(video, i);
    frame_temp = im2double(rgb2gray(frame_temp));
    frame_dif = frame_temp - b_pre; 
    b_temp = b_pre;
    b_temp(abs(frame_dif) < 3*std_mat) = b_pre(abs(frame_dif) < 3*std_mat)*p + frame_temp(abs(frame_dif) < 3*std_mat)*(1-p); % 当前背景
    % imshow(b_temp)
    frame_image = frame_temp - b_temp;
    % imshow(frame_image,[]);
    % imshow(im2bw(abs(frame_image)))
%     imshow(abs(frame_image),[])
    frame_image = im2bw(abs(frame_image),0.2); %二值化显示前景
    B = [0 1 0; 1 1 1; 0 1 0];
    frame_image = imdilate(frame_image, B); %膨胀操作
%     frame_image = imdilate(frame_image, B);
    imshow(frame_image)
    b_pre = b_temp; % 更新背景   
end
    
    
