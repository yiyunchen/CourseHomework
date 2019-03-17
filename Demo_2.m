%% ���˹������ģ
clear;
video = VideoReader('768x576.avi');
nFrames = video.NumberOfFrames;
fr = read(video, 1);
% ��ȡ��ͼ����Ϊ���� 
fr_bw1 = rgb2gray(fr);   
% ������ת��Ϊ�Ҷ�ͼ�� 
fr_size = size(fr);             
%ȡ֡��С 
width = fr_size(2); 
height = fr_size(1); %��ȡԭͼ��ĳߴ�
fg = zeros(height, width);   
bg_bw = fr_bw1;
fr_bw1 = double(fr_bw1);

C = 3;           % ��ɻ�ϸ�˹�ĵ���˹��Ŀ (һ��3-5) 
D = 2.5;         % ��ֵ��һ��2.5����׼� 
alpha = 0.01;    % learning rate ѧϰ�ʾ��������ٶ�(between 0 and 1) 
% thresh = 0.25;   % foreground threshold ǰ����ֵ
sd_init = 36;    % initial standard deviation ��ʼ����׼��(for new components)
w = ones(height,width,C);          % initialize weights array ��ʼ��Ȩֵ���� 
w(:,:,1) = 0.9;
w(:,:,2:C) = 0.1;                  % ��һ����˹�ֲ��ĳ�ʼȨ��Ϊ1������ֲ���Ȩ��Ϊ0
avg = zeros(height,width,C);        % pixel means ���ؾ�ֵ 
avg(:,:,1) = fr_bw1;                % ��һ����˹�ֲ��ĳ�ʼ��ֵΪ�ο�֡��ֵ������ֲ��ľ�ֵΪ0s
sd = sd_init*ones(height,width,C);   % pixel standard deviations ���ر�׼�� 
match = zeros(height, width,C);    % ƥ��Ĵ�������ʼֵ����Ϊ0
% matchcnt(:,:,1) = ones(height, width);
u_diff = zeros(height,width,C);      % difference of each pixel from mean ͼƬ���˹��ֵ�Ĳ� 

for num = 2:nFrames
    fr = read(video, num);
    fr_bw = rgb2gray(fr);  	   % convert frame to grayscale ת��Ϊ�Ҷ�ͼ��
    fr_bw = double(fr_bw);     % ���Ҷ�ͼֵ����Ϊ˫����
    %���������ͼƬ�������˹��ֵ�Ĳ�
    for m=1:C  
        u_diff(:,:,m) = abs(fr_bw - avg(:,:,m));
        match(:,:,m) = u_diff(:,:,m) < 2.5*sd(:,:,m);
        w_temp = w(:,:,m);
        w_temp(match(:,:,m) == 1) = (1 - alpha) * w_temp(match(:,:,m) == 1) + alpha;
        w_temp(match(:,:,m) == 0) = w_temp(match(:,:,m) == 0)*(1 - alpha);
        w(:,:,m) = w_temp;
        p = alpha./w_temp;
        avg_temp = avg(:,:,m);
        avg_temp(match(:,:,m) == 1) = (1 - p(match(:,:,m) == 1)) .* avg_temp(match(:,:,m) == 1) + p(match(:,:,m) == 1) .* fr_bw(match(:,:,m) == 1);
        avg(:,:,m) = avg_temp;
        sd_temp = sd(:,:,m);      
        sd_temp(match(:,:,m) == 1) = sqrt((1 - p(match(:,:,m) == 1)) .* sd_temp(match(:,:,m) == 1).^2 + p(match(:,:,m) == 1) .* (fr_bw(match(:,:,m) == 1) - avg(match(:,:,m) == 1)).^2);
        sd(:,:,m) = sd_temp;     
    end      
    w_sum = sum(w, 3);
    for m=1:C 
        w(:,:,m) = w(:,:,m)./w_sum;
    end
    
    rank = w./sd; % ����
    [rank index] = sort(rank,3,'descend');
    for j = 1:width
        for i = 1:height
            temp = [avg(i,j,:); sd(i,j,:); w(i,j,:)];
            avg(i,j,:) = temp(1,index(i,j,:));
            sd(i,j,:) = temp(2,index(i,j,:));
            w(i,j,:) = temp(3,index(i,j,:));
        end
    end
    
    match_sum = sum(match,3);
    [m_vec, n_vec] = find(match_sum == 0);
    avg(m_vec,n_vec,3) = fr_bw(m_vec,n_vec);
    sd(m_vec,n_vec,3) = sd_init;
    
    w_cumsum = cumsum(w,3);
    image = 255*ones(height,width); % ǰ��
    for j = 1:width
        for i = 1:height
            w_cumsum_temp = w_cumsum(i,j,:);
            index = find(w_cumsum_temp > 0.7);
            flag = 0;
            for k = 1:index(1) % �����ж�
                if u_diff(i,j,k) < 2.5*sd(i,j,k)
                    flag = 1;
                end
            end
            if flag == 1
                image(i,j) = 0;
            end
        end
    end
    B = [0 1 0; 1 1 1; 0 1 0];
    image = imdilate(image, B); %���Ͳ���
    imshow(image,[])
    match = zeros(height, width, C);    
%     for i = 1:length(m_vec)
%         for j = 1:length(n_vec)
% %             [min_w, min_w_index] = min(w(m_vec(i),n_vec(j),:));
%             avg(m_vec(i),n_vec(j),min_w_index) = fr_bw(m_vec(i),n_vec(j));
%             sd(m_vec(i),n_vec(j),min_w_index) = sd_init;
%         end
%     end   
%     bg_bw = 0.9 * bg_bw + 0.1 * sum(avg.*w,3);
end

				

