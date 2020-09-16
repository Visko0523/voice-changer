function varargout = GUI(varargin)
% GUI MATLAB code for GUI.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI

% Last Modified by GUIDE v2.5 09-Sep-2020 21:11:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before GUI is made visible.
function GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI (see VARARGIN)

handles.Fs=8000; % sampling frequency

% Choose default command line output for GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in rec.
function rec_Callback(hObject, eventdata, handles)
% hObject    handle to rec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    global r
    Fs=handles.Fs;
    if get(hObject,'value')
        % start recodring
        r = audiorecorder(Fs, 16, 1);
        record(r); % record data
    else
        % stop recording
        stop(r);
        s = getaudiodata(r); % get data
        s=s-mean(s);
        handles.recs=s';
    end
    guidata(hObject,handles);

% Hint: get(hObject,'Value') returns toggle state of rec


% --- Executes on button press in yuansheng.
function yuansheng_Callback(hObject, eventdata, handles)
% hObject    handle to yuansheng (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    recs=handles.recs;
    Fs=handles.Fs;
    sound(recs,Fs)%����ԭ��
    
    N=length(recs);%��ȡ��Ƶ����
    F=abs(fft(recs))/N*2;%����Ƶ�׷���
    f=(0:N-1)*Fs/N;%Ƶ�׺����꣨Hz��
    t=(0:N-1)/Fs;%ʱ������꣨s��
    M=find(floor(f/3400),1,'first');%��ѡ3400Hz��������Χ�����µ�Ƶ��
    
    plot(handles.axes1,t,recs)
    title(handles.axes1,'ʱ����');
    xlabel(handles.axes1,'t/s');ylabel(handles.axes1,'��ֵ');

    plot(handles.axes2,f(1:M),F(1:M))
    title(handles.axes2,'��Ƶ����');
    xlabel(handles.axes2,'f/Hz');ylabel(handles.axes2,'��ֵ');
 

% --- Executes on button press in nvbiannan.
function nvbiannan_Callback(hObject, eventdata, handles)
% hObject    handle to nvbiannan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    x=handles.recs;
    Fs=handles.Fs;
    M=length(x);
    xaa=fft(x);
    pa=[xaa(1:floor(0.8*M)),xaa(floor(0.5*M):M)];
    Y1=5*real(ifft(pa));
    sound(Y1,Fs);
    N=length(Y1);
    F=abs(fft(Y1))/N*2;%����Ƶ�׷���
    f=(0:N-1)*Fs/N;%Ƶ�׺����꣨Hz��
    t=(0:N-1)/Fs;%ʱ������꣨s��
    M=find(floor(f/3400),1,'first');%��ѡ3400Hz��������Χ�����µ�Ƶ��
    plot(handles.axes1,t,Y1)
    title(handles.axes1,'ʱ����');
    xlabel(handles.axes1,'t/s');ylabel(handles.axes1,'��ֵ');

    plot(handles.axes2,f(1:M),F(1:M))
    title(handles.axes2,'��Ƶ����');
    xlabel(handles.axes2,'f/Hz');ylabel(handles.axes2,'��ֵ');



% --- Executes on button press in nanbiannv.
function nanbiannv_Callback(hObject, eventdata, handles)
% hObject    handle to nanbiannv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    s=handles.recs;
    Fs=handles.Fs;
    FL = 80 ; % ֡��
    WL = 240 ; % ����
    P = 10 ; %Ԥ��ϵ������
    L = length(s); % ������������
    FN = floor(L/FL)-2; % ����֡����floor��������� 
    % Ԥ����ؽ��˲���
    exc = zeros(L,1); % �����źţ�double�������L��1��
    zi_pre = zeros(P,1); % Ԥ���˲���״̬
    s_rec = zeros(L,1); % �ؽ�����
    zi_rec = zeros(P,1); % ����˲���
    exc_syn_t = zeros(L,1); % �ϳɵļ����źţ�����һ��L��1�е�0����s_syn_t = zeros(L,1); % �ϳ�����
    last_syn_t = 0; % �洢��һ���ε����һ��������±�
    zi_syn_t = zeros(P,1); % �ϳ��˲���
    hw = hamming(WL); %������
    %�˲���
    % ���δ���ÿ֡����
    for n = 3:FN %�ӵ����������鿪ʼ
        % ����Ԥ��ϵ��
        s_w = s(n*FL-WL+1:n*FL)'.*hw; %��������Ȩ
        [A,E]=lpc(s_w,P); %����Ԥ�����Ԥ��ϵ����LPCϵ����
        % A��Ԥ��ϵ����E�ᱻ��������ϳɼ���������
        s_f=s((n-1)*FL+1:n*FL); % ��֡����
        %����filter�����ؽ�����
        [exc1,zi_pre] = filter(A,1,s_f,zi_pre);
        exc((n-1)*FL+1:n*FL) = exc1; %���㼤��
        %����filter�����ؽ�����
        [s_rec1,zi_rec] = filter(1,A,exc1,zi_rec);
        s_rec((n-1)*FL+1:n*FL) = s_rec1; %�ؽ�����
        % ����ֻ�еõ�exc��ſ���
        s_Pitch = exc(n*FL-222:n*FL);

        PT = findpitch(s_Pitch); %�����������pt
        G = sqrt(E*PT); %����ϳɼ���������G

        PT1 =floor(PT/3);
        poles = roots(A);
        deltaOMG =300*2*pi/8000;
        for p=1:10
            if imag(poles(p))>0 ;
            poles(p) = poles(p)*exp(1i*deltaOMG);
            elseif imag(poles(p))<0 ;
            poles(p) = poles(p)*exp(-1i*deltaOMG);
            end
        end
        A1=poly(poles);
        tempn_syn_t = 1:n*FL-last_syn_t;
        exc_syn1_t = zeros(length(tempn_syn_t),1);
        exc_syn1_t(mod(tempn_syn_t,PT1)==0) = G;
        exc_syn1_t = exc_syn1_t((n-1)*FL-last_syn_t+1:n*FL-last_syn_t); 
        [s_syn1_t,zi_syn_t] = filter(1,A1,exc_syn1_t,zi_syn_t);
        exc_syn_t((n-1)*FL+1:n*FL) = exc_syn1_t;
        s_syn_t((n-1)*FL+1:n*FL) = s_syn1_t;
        last_syn_t = last_syn_t+PT1*floor((n*FL-last_syn_t)/PT1);
    end
    y=3*s_syn_t;
    sound(y,Fs);
    N=length(y);
    F=abs(fft(y))/N*2;%����Ƶ�׷���
    f=(0:N-1)*Fs/N;%Ƶ�׺����꣨Hz��
    t=(0:N-1)/Fs;%ʱ������꣨s��
    M=find(floor(f/3400),1,'first');%��ѡ3400Hz��������Χ�����µ�Ƶ��
    plot(handles.axes1,t,y)
    title(handles.axes1,'ʱ����');
    xlabel(handles.axes1,'t/s');ylabel(handles.axes1,'��ֵ');
    plot(handles.axes2,f(1:M),F(1:M))
    title(handles.axes2,'��Ƶ����');
    xlabel(handles.axes2,'f/Hz');ylabel(handles.axes2,'��ֵ');

% --- Executes on button press in tongsheng.
function tongsheng_Callback(hObject, eventdata, handles)
% hObject    handle to tongsheng (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    s=handles.recs;
    Fs=handles.Fs;
    FL = 80 ; % ֡��
    WL = 240 ; % ����
    P = 10 ; %Ԥ��ϵ������
    L = length(s); % ������������
    FN = floor(L/FL)-2; % ����֡����floor��������� 
    % Ԥ����ؽ��˲���
    exc = zeros(L,1); % �����źţ�double�������L��1��
    zi_pre = zeros(P,1); % Ԥ���˲���״̬
    s_rec = zeros(L,1); % �ؽ�����
    zi_rec = zeros(P,1); % ����˲���
    exc_syn_t = zeros(L,1); % �ϳɵļ����źţ�����һ��L��1�е�0����s_syn_t = zeros(L,1); % �ϳ�����
    last_syn_t = 0; % �洢��һ���ε����һ��������±�
    zi_syn_t = zeros(P,1); % �ϳ��˲���
    hw = hamming(WL); %������
    %�˲���
    % ���δ���ÿ֡����
    for n = 3:FN %�ӵ����������鿪ʼ
        % ����Ԥ��ϵ��
        s_w = s(n*FL-WL+1:n*FL)'.*hw; %��������Ȩ
        [A,E]=lpc(s_w,P); %����Ԥ�����Ԥ��ϵ����LPCϵ����
        % A��Ԥ��ϵ����E�ᱻ��������ϳɼ���������
        s_f=s((n-1)*FL+1:n*FL); % ��֡����
        %����filter�����ؽ�����
        [exc1,zi_pre] = filter(A,1,s_f,zi_pre);
        exc((n-1)*FL+1:n*FL) = exc1; %���㼤��
        %����filter�����ؽ�����
        [s_rec1,zi_rec] = filter(1,A,exc1,zi_rec);
        s_rec((n-1)*FL+1:n*FL) = s_rec1; %�ؽ�����
        % ����ֻ�еõ�exc��ſ���
        s_Pitch = exc(n*FL-222:n*FL);

        PT = findpitch(s_Pitch); %�����������pt
        G = sqrt(E*PT); %����ϳɼ���������G

        PT1 =floor(PT/2);
        poles = roots(A);
        deltaOMG =300*2*pi/8000;
        for p=1:10
            if imag(poles(p))>0 ;
            poles(p) = poles(p)*exp(1i*deltaOMG);
            elseif imag(poles(p))<0 ;
            poles(p) = poles(p)*exp(-1i*deltaOMG);
            end
        end
        A1=poly(poles);
        tempn_syn_t = 1:n*FL-last_syn_t;
        exc_syn1_t = zeros(length(tempn_syn_t),1);
        exc_syn1_t(mod(tempn_syn_t,PT1)==0) = G;
        exc_syn1_t = exc_syn1_t((n-1)*FL-last_syn_t+1:n*FL-last_syn_t); 
        [s_syn1_t,zi_syn_t] = filter(1,A1,exc_syn1_t,zi_syn_t);
        exc_syn_t((n-1)*FL+1:n*FL) = exc_syn1_t;
        s_syn_t((n-1)*FL+1:n*FL) = s_syn1_t;
        last_syn_t = last_syn_t+PT1*floor((n*FL-last_syn_t)/PT1);
    end
    y=3*s_syn_t;
    sound(y,Fs);
    N=length(y);
    F=abs(fft(y))/N*2;%����Ƶ�׷���
    f=(0:N-1)*Fs/N;%Ƶ�׺����꣨Hz��
    t=(0:N-1)/Fs;%ʱ������꣨s��
    M=find(floor(f/3400),1,'first');%��ѡ3400Hz��������Χ�����µ�Ƶ��
    plot(handles.axes1,t,y)
    title(handles.axes1,'ʱ����');
    xlabel(handles.axes1,'t/s');ylabel(handles.axes1,'��ֵ');
    plot(handles.axes2,f(1:M),F(1:M))
    title(handles.axes2,'��Ƶ����');
    xlabel(handles.axes2,'f/Hz');ylabel(handles.axes2,'��ֵ');


% --- Executes on button press in laoren.
function laoren_Callback(hObject, eventdata, handles)
% hObject    handle to laoren (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    x=handles.recs;
    Fs=handles.Fs;
    M=length(x);

    %��ͨ�˲��������
    fp1=2500;
    fs1=2600;
    %�趨��ͨ�˲���ͨ����ֹƵ�ʺ������ֹƵ��
    wp1=2*fp1/Fs;
    ws1=2*fs1/Fs;
    rp=1;
    as=100;
    [N1,wp1]=ellipord(wp1,ws1,rp,as); %�����ͨ�˲���������ͨ���߽�Ƶ��
    [B,A]=ellip(N1,rp,as,wp1); %�����ͨ�˲���ϵͳ����ϵ��
    y1=filter(B,A,x);

    xaa=fft(y1);

    pa=[xaa(1:floor(0.9*M)),xaa(floor(0.3*M):M)];
    Y1=5*real(ifft(pa));
    sound(Y1,Fs);
    N=length(Y1);
    F=abs(fft(Y1))/N*2;%����Ƶ�׷���
    f=(0:N-1)*Fs/N;%Ƶ�׺����꣨Hz��
    t=(0:N-1)/Fs;%ʱ������꣨s��
    M=find(floor(f/3400),1,'first');%��ѡ3400Hz��������Χ�����µ�Ƶ��
    plot(handles.axes1,t,Y1)
    title(handles.axes1,'ʱ����');
    xlabel(handles.axes1,'t/s');ylabel(handles.axes1,'��ֵ');

    plot(handles.axes2,f(1:M),F(1:M))
    title(handles.axes2,'��Ƶ����');
    xlabel(handles.axes2,'f/Hz');ylabel(handles.axes2,'��ֵ');


% --- Executes on button press in xiaohuang.
function xiaohuang_Callback(hObject, eventdata, handles)
% hObject    handle to xiaohuang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    recs=handles.recs;
    Fs=1.8*handles.Fs;
    sound(recs,Fs)
    
    N=length(recs);%��ȡ��Ƶ����
    F=abs(fft(recs))/N*2;%����Ƶ�׷���
    f=(0:N-1)*Fs/N;%Ƶ�׺����꣨Hz��
    t=(0:N-1)/Fs;%ʱ������꣨s��
    M=find(floor(f/3400),1,'first');%��ѡ3400Hz��������Χ�����µ�Ƶ��
    
    plot(handles.axes1,t,recs)
    title(handles.axes1,'ʱ����');
    xlabel(handles.axes1,'t/s');ylabel(handles.axes1,'��ֵ');

    plot(handles.axes2,f(1:M),F(1:M))
    title(handles.axes2,'��Ƶ����');
    xlabel(handles.axes2,'f/Hz');ylabel(handles.axes2,'��ֵ');


% --- Executes on button press in shulan.
function shulan_Callback(hObject, eventdata, handles)
% hObject    handle to shulan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    recs=2*handles.recs;
    Fs=0.5*handles.Fs;
    sound(recs,Fs)
    
    N=length(recs);%��ȡ��Ƶ����
    F=abs(fft(recs))/N*2;%����Ƶ�׷���
    f=(0:N-1)*Fs/N;%Ƶ�׺����꣨Hz��
    t=(0:N-1)/Fs;%ʱ������꣨s��
    M=find(floor(f/3400),1,'first');%��ѡ3400Hz��������Χ�����µ�Ƶ��
    
    plot(handles.axes1,t,recs)
    title(handles.axes1,'ʱ����');
    xlabel(handles.axes1,'t/s');ylabel(handles.axes1,'��ֵ');

    plot(handles.axes2,f(1:M),F(1:M))
    title(handles.axes2,'��Ƶ����');
    xlabel(handles.axes2,'f/Hz');ylabel(handles.axes2,'��ֵ');
