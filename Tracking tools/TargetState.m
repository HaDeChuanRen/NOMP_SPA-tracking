classdef TargetState
    %   last date: 2023/7/20
    %
    %   2023/7/20 update detials:
    %   1. (Todo) The probability of detection (P_det) is substituted for appear
    %    state variable (appear_state) to repensent the detection state of
    %    targets.
    %   2. (Todo) The probability of detection (P_det) increases when the
    %     poseterior state is estimated, decreases when the target disappearing.
    %   3. (Todo) rewrite and add the comments of the classdef file.
    %
    %   TargetState: Summary of this class goes here
    %   Detailed explanation goes here
    %
    %   Properties
    %   statevec: state vector [Px, Vx, Py, Vy]^{\rm T} denotes x position, y
    %   position, x velocity,  y velocity of target, respectively
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
            obj.P_det = 0.5;
            obj.PD_his = obj.P_det;
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

        function obj = state_estimate(obj, Rmat_noise, Hmat, ZPost)
            % obj.statevec = obj.statevec + Kmat_gain * (ZPost' - Hmat * ...
            %     obj.statevec);
            % obj.Pmat = (eye(obj.state_dim) - Kmat_gain * Hmat) * obj.Pmat;
            Smat_t = Hmat * obj.Pmat * Hmat' + Rmat_noise;
            I_dim = eye(obj.state_dim);
            Pri_mat = obj.Pmat;
            obj.Pmat = I_dim / (Hmat' / Smat_t * Hmat + I_dim / Pri_mat);
            obj.statevec = obj.Pmat * (Hmat' / Smat_t * ZPost' + ...
                Pri_mat \ obj.statevec);
            obj.P_det = obj.P_det * 1.25;
            if obj.P_det > 1
                obj.P_det = 1;
            end

            obj.PD_his = [obj.PD_his; obj.P_det];
            state_history = obj.state_his;
            obj.last_time = obj.last_time + 1;
            obj.state_his = zeros(obj.state_dim, obj.last_time);
            obj.state_his(:, 1 : end - 1) = state_history;
            obj.state_his(:, end) = obj.statevec;
        end

        function obj = disappearing(obj)
            obj.P_det = obj.P_det / 1.25;
            if obj.P_det > 0.1
                obj.last_time = obj.last_time + 1;
                % update history
                state_history = obj.state_his;
                obj.state_his = zeros(obj.state_dim, obj.last_time);
                obj.state_his(:, 1 : end - 1) = state_history;
                obj.state_his(:, end) = obj.statevec;
                obj.PD_his = [obj.PD_his; obj.P_det];
            end
        end
    end
end