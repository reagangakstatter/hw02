% Author: Reagan Gakstatter / reg0052@auburn.edu
% Date: 2024-10-01
% Assignment Name: hw02

classdef hw02
    methods (Static)

        function [c, n] = p1(f, a, b, epsilon, name, f_prime)
            % p1: Implement numerical methods to find the root of a function
            % :param f: function handle
            % :param a: real number, the left end of the interval
            % :param b: real number, the right end of the interval
            % :param epsilon: real number, the function tolerance
            % :param name: string, the name of the method
            % :param f_prime: function handle, the derivative of the function, only needed for Newton's method
            %
            % :return: 
            % c: real number, the root of the function
            % n: integer, the number of iterations

            % Bisection method
            if strcmp(name, 'bisection')
                n = 0;
                c = (a + b) / 2;
                while abs(f(c)) > epsilon
                    if f(a) * f(c) < 0
                        b = c;
                    else
                        a = c;
                    end
                    c = (a + b) / 2;
                    n = n + 1;
                end

            % Secant method
            elseif strcmp(name, 'secant')
                n = 0;
                c = a;
                c_prev = b;
                while abs(f(c)) > epsilon
                    temp = c;
                    c = c - f(c) * (c - c_prev) / (f(c) - f(c_prev));
                    c_prev = temp;
                    n = n + 1;
                end

            % Newton's method
            elseif strcmp(name, 'newton')
                if nargin < 6
                    error('f_prime (the derivative of f) is required for Newton''s method.');
                end
                n = 0;
                c = a; % Initial guess
                while abs(f(c)) > epsilon
                    c = c - f(c) / f_prime(c);
                    n = n + 1;
                end

            % False position (Regula Falsi) method
            elseif strcmp(name, 'regula_falsi')
                n = 0;
                while abs(b - a) > epsilon
                    c = a - (f(a) * (b - a)) / (f(b) - f(a));
                    if abs(f(c)) < epsilon
                        break;
                    end
                    if f(a) * f(c) < 0
                        b = c;
                    else
                        a = c;
                    end
                    n = n + 1;
                end

            % Steffensen's method
            elseif strcmp(name, 'steffensen')
                n = 0;
                c = a; % Initial guess
                while abs(f(c)) > epsilon
                    g = f(c);
                    c_next = c - (g^2 / (f(c + g) - g));
                    c = c_next;
                    n = n + 1;
                end
            else
                error('Unknown method name.');
            end
        end

        function p2()
            % Summarize the iteration number for each method in the table
            %     |name          | iter | 
            %     |--------------|------|
            %     |bisection     |32    |
            %     |secant        |8     |
            %     |newton        |5     |
            %     |regula_falsi  |59    |
            %     |steffensen    |13    |
            f = @(x) 3*x^3 - 2*x^2 - 4;
            f_prime = @(x) 9*x^2 - 4*x;
            a = 1; b = 3;
            epsilon = 1e-9;

            methods = {'bisection', 'secant', 'newton', 'regula_falsi', 'steffensen'};
            fprintf('|%-15s|%-5s|\n', 'Method', 'Iter');
            fprintf('|---------------|-----|\n');
            for i = 1:length(methods)
                [~, n] = hw02.p1(f, a, b, epsilon, methods{i}, f_prime);
                fprintf('|%-15s|%-5d|\n', methods{i}, n);
            end
        end

    end
end
