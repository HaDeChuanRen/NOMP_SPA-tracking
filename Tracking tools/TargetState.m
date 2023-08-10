classdef TargetState
    %   last date: 2023/7/24
    %   2023/7/24 update detials:
    %   1. (Todo) rewrite and add the comments of the classdef file.
    %   2. (Todo) add default parameters for the functions.
    %   3. (Todo) estimate the posterier target location by weighting the
    %       measurements with the correlation probability.
    %   4. (Todo) use the appear state flag to denote whether the target appear %       or not.
    %
    %   Properties
    %   statevec: state vector [Px, Vx, Py, Vy]^{\rm T} denotes x position, y
    %       position, x velocity,  y velocity of target, respectively
    %   state_dim: the dimension of statevec, 4 in genenral
    %   n_label: the label of target
    %   Pmat: the covariance matrix of Px, Py, Vx, Vy
    %   Qmat: transition covariance ???
    %   Amat: transition matrix
    %   cu_moment: current moment
    %   T_interval: time intervel (sampling interval time)
    %   last_time: the target last time
    %   state_his: the history x position, y position, x velocity, y velocity
    %   P_det: the detection probability of target.
    %   PD_his: the history of detection probability.
    %
    %   Methods

    properties
        statevec
        state_dim
        n_label
        Pmat
        Qmat
        Amat
        cu_moment
        T_interval
        last_time
        state_his
        disapp_his
        disapp_time
        app_state
        P_det
        PD_his
    end

    methods
        function obj = TargetState(instate, inlabel, inPmat, inQmat, ...
            inAmat, in_moment, in_interval)
            %   Construct an instance of this class
            %   The input parameter:
            %   instate: the input of original target state
            %   inlabel: the input of target label
            %   inPmat: the input of the covariance matrix of Px, Py, Vx, Vy
            %   inQmat: the input of transition covariance
            %   inAmat: the input of transition matrix
            %   in_moment: the input of current moment
            %   in_interval: the input of observation interval

            obj.statevec = instate;
            obj.state_dim = length(instate);
            obj.n_label = inlabel;
            obj.Pmat = inPmat;
            obj.Qmat = inQmat;
            obj.Amat = inAmat;
            obj.cu_moment = in_moment;
            obj.T_interval = in_interval;
            obj.last_time = 1;
            obj.state_his = obj.statevec;
            obj.P_det = 0.8;
            obj.PD_his = obj.P_det;
            obj.disapp_his = [];
            obj.disapp_time = 0;
            obj.app_state = 1;
        end

        function [state_new, Pmat_pri, obj] = state_transform(obj, inAmat, in_T)
            %   Summary of this method goes here
            %   Detailed explanation goes here
            if ~exist('guard_band','var'), inAmat = obj.Amat;
            elseif isempty(inAmat), inAmat = obj.Amat; end
            if ~exist('guard_band','var'), in_T = obj.T_interval;
            elseif isempty(inAmat), in_T = obj.T_interval; end
            obj.Amat = inAmat;
            obj.T_interval = in_T;

            state_new = obj.Amat * obj.statevec;
            Pmat_pri = obj.Amat * obj.Pmat * obj.Amat' + obj.Qmat;
            obj.cu_moment = obj.cu_moment + obj.T_interval;
            % update state
            obj.statevec = state_new;
            obj.Pmat = Pmat_pri;
            % update history
            % obj.last_time = obj.last_time + 1;
            % state_history = obj.state_his;
            % obj.state_his = zeros(obj.state_dim, obj.last_time);
            % obj.state_his(:, 1 : end - 1) = state_history;
            % obj.state_his(:, end) = obj.statevec;

        end

        function obj = est_hard(obj, Rmat_noise, Hmat, ZPost)
            % obj.statevec = obj.statevec + Kmat_gain * (ZPost' - Hmat * ...
            %     obj.statevec);
            % obj.Pmat = (eye(obj.state_dim) - Kmat_gain * Hmat) * obj.Pmat;
            % Smat_t = Hmat * obj.Pmat * Hmat' + Rmat_noise;
            I_dim = eye(obj.state_dim);
            Pri_mat = obj.Pmat;
            state_pri = obj.statevec;
            Kmat_gain = Pri_mat * Hmat' / (Hmat * Pri_mat * ...
                Hmat' + Rmat_noise);

            Post_mat = (I_dim - Kmat_gain * Hmat) * Pri_mat;
            state_post = state_pri + Kmat_gain * (ZPost' - Hmat * state_pri);
            % obj.Pmat = I_dim / (Hmat' / Smat_t * Hmat + I_dim / Pri_mat);
            % obj.statevec = obj.Pmat * (Hmat' / Smat_t * ZPost' + ...
            %     Pri_mat \ obj.statevec);

            % update state and history
            obj.Pmat = Post_mat;
            obj.statevec = state_post;
            obj.P_det = obj.P_det * 1.1;
            if obj.P_det > 0.999
                obj.P_det = 0.999;
            end
            obj.PD_his = [obj.PD_his; obj.P_det];

            if obj.disapp_time == 0
                state_history = obj.state_his;
                obj.last_time = obj.last_time + 1;
                obj.state_his = zeros(obj.state_dim, obj.last_time);
                obj.state_his(:, 1 : end - 1) = state_history;
                obj.state_his(:, end) = obj.statevec;
            else
                state_history = obj.state_his;
                old_last = obj.last_time;
                obj.last_time = obj.last_time + obj.disapp_time + 1;
                obj.state_his = zeros(obj.state_dim, obj.last_time);
                obj.state_his(:, 1 : old_last) = state_history;
                obj.state_his(:, (old_last + 1) : (end - 1)) = obj.disapp_his;
                obj.state_his(:, end) = obj.statevec;

                obj.disapp_his = [];
                obj.disapp_time = 0;
            end
        end

        function obj = est_soft(obj, Rmat_noise, Hmat, Zset, pro_vec)
            I_dim = eye(obj.state_dim);
            Pri_mat = obj.Pmat;
            state_pri = obj.statevec;
            Kmat_gain = Pri_mat * Hmat' / (Hmat * Pri_mat * ...
                Hmat' + Rmat_noise);

            ZPri_kt = Hmat * state_pri;
            Z_ave = [Zset; ZPri_kt']' * pro_vec;
            Post_mat = (I_dim - Kmat_gain * Hmat) * Pri_mat;
            state_post = state_pri + Kmat_gain * (Z_ave - Hmat * state_pri);
            obj.Pmat = Post_mat;
            obj.statevec = state_post;

            if pro_vec(end) < 0.5
                obj.P_det = obj.P_det * 1.1;
                if obj.P_det > 0.999
                    obj.P_det = 0.999;
                end
                obj.PD_his = [obj.PD_his; obj.P_det];

                if obj.disapp_time == 0
                    state_history = obj.state_his;
                    obj.last_time = obj.last_time + 1;
                    obj.state_his = zeros(obj.state_dim, obj.last_time);
                    obj.state_his(:, 1 : end - 1) = state_history;
                    obj.state_his(:, end) = obj.statevec;
                else
                    state_history = obj.state_his;
                    old_last = obj.last_time;
                    obj.last_time = obj.last_time + obj.disapp_time + 1;
                    obj.state_his = zeros(obj.state_dim, obj.last_time);
                    obj.state_his(:, 1 : old_last) = state_history;
                    obj.state_his(:, (old_last + 1) : (end - 1)) = obj.disapp_his;
                    obj.state_his(:, end) = obj.statevec;

                    obj.disapp_his = [];
                    obj.disapp_time = 0;
                end
            else
                obj = obj.disappearing();
            end
        end

        function obj = disappearing(obj)
            obj.P_det = obj.P_det * 0.9;
            if obj.P_det > 0.5
                % obj.last_time = obj.last_time + 1;
                % % update history
                % state_history = obj.state_his;
                % obj.state_his = zeros(obj.state_dim, obj.last_time);
                % obj.state_his(:, 1 : end - 1) = state_history;
                % obj.state_his(:, end) = obj.statevec;
                obj.disapp_time = obj.disapp_time + 1;
                dis_history = obj.disapp_his;
                obj.disapp_his = zeros(obj.state_dim, obj.disapp_time);
                obj.disapp_his(:, 1 : end - 1) = dis_history;
                obj.disapp_his(:, end) = obj.statevec;
                obj.PD_his = [obj.PD_his; obj.P_det];
            else
                obj.app_state = 0;
            end
        end
    end
end

%   2023/7/22 update detials:
%   1. (Todo) rewrite and add the comments of the classdef file.
%   2. (Todo) add default parameters for the functions.
%   3. revise the Kalman filter
%   4. seperate the histoty of the stage when the target is
%       disappearing. If there is new information, the hisroty of
%       disappearing will be written into the total history.
%
%   2023/7/20 update detials:
%   1. The probability of detection (P_det) is substituted for appear
%    state variable (appear_state) to repensent the detection state of
%    targets.
%   2. The probability of detection (P_det) increases when the
%     poseterior state is estimated, decreases when the target disappearing.
%
%   TargetState: Summary of this class goes here
%   Detailed explanation goes here