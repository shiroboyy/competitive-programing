#include <iostream>
#include <string>
#include <cmath>
using namespace std;

float round_to_decimal(float var, int decimal_places) {
    float factor = pow(10.0f, decimal_places);
    return round(var * factor) / factor;
}

// 6.1
long long compute_factorial(int n);
long long compute_combination(int n, int k);

// 6.2
double exponentiation(double x, int n);
double compute_equation(double x, int n);

// 6.3
double power_of_e(double x, int n = 20);

// 6.4
struct QuadraticSolution {
    bool has_solution;
    double x1, x2;
};
QuadraticSolution solve_quadratic(double a, double b, double c);

// 6.5
struct EquationSolution {
    int status;
    double x, y;
};
EquationSolution solve_system_equation(double a, double b, double c,
                                       double d, double e, double f);

// 6.6
string format_read_number(string number_str) {
    string result = "";
    for (int i = 0; i < number_str.length(); i++) {
        if (isspace(number_str[i]) && isspace(result.back())){
            continue;
        } else {
            result += number_str[i];
        }
    }
    if (isspace(result.back())) {
        result.pop_back();
    }
    return result;
}
string read_number(int n);

// 6.8
float convert_temperature(float temperature, bool f2c = true);

// 6.9
struct Point {
    float x;
    float y;
};
float linear_regression(Point data_points[], int num_data_points, float new_data,
                        float learning_rate = 1e-3, int num_iterations = 5000);


int main() {
    cout << compute_combination(9, 5) << "\n";
    cout << compute_equation(0.2, 3) << "\n";
    cout << power_of_e(3) << "\n";
    QuadraticSolution sol = solve_quadratic(4, -2, -6); 
    cout << sol.x1 << " " << sol.x2 << "\n";
    EquationSolution sol_1 = solve_system_equation(-2, 1, 5, 1, 3, 1);
    cout << sol_1.x << " " << sol_1.y << "\n"; 
    cout << read_number(900000000) << endl;
    cout << read_number(900000021) << endl;
    cout << read_number(900100021) << endl;
    cout << read_number(900001021) << endl;
    cout << read_number(900021001) << endl;
    cout << read_number(900000001) << endl;
    cout << read_number(900000100) << endl;
    cout << convert_temperature(33, false) << "\n";

    Point points[5] = {{1, 2}, {3, 6}, {2, 4}, {5, 10}, {11, 22}};
    cout << linear_regression(points, 5, 6) << "\n";
    return 0;
}


long long compute_factorial(int n) {
    long long ans = 1;
    for (int i = 1; i <= n; i++)
        ans *= i;
    return ans;
}

long long compute_combination(int n, int k) {
    return compute_factorial(n) /
           (compute_factorial(k) * compute_factorial(n - k));
}

double exponentiation(double x, int n) {
    double ans = 1;
    for (int i = 1; i <= n; i++)
        ans *= x;
    return ans;
}

double compute_equation(double x, int n) {
    double sum = exponentiation(1.5, 8);
    for (int i = 1; i <= n; i++) {
        double term = exponentiation(x + i, i) / (i * i);
        sum += (i % 2 == 0 ? term : -term);
    }
    return round_to_decimal(sum, 3);
}

double power_of_e(double x, int n) {
    double result = 1.0;
    for (int i = 1; i <= n; i++)
        result += exponentiation(x, i) / compute_factorial(i);
    return round_to_decimal(result, 5);
}

QuadraticSolution solve_quadratic(double a, double b, double c) {
    QuadraticSolution res = {false, 0, 0};
    double delta = b * b - 4 * a * c;
    if (delta < 0) return res;
    if (delta == 0) {
        double x = round_to_decimal(-b / (2 * a), 3);
        res = {true, x, x};
    } else {
        double x1 = round_to_decimal((-b + sqrt(delta)) / (2 * a), 3);
        double x2 = round_to_decimal((-b - sqrt(delta)) / (2 * a), 3);
        if (x1 > x2) swap(x1, x2);
        res = {true, x1, x2};
    }
    return res;
}

EquationSolution solve_system_equation(double a, double b, double c,
                                       double d, double e, double f) {
    EquationSolution res;
    double detD = a * e - d * b;
    double detDx = c * e - f * b;
    double detDy = a * f - d * c;
    if (detD != 0) {
        res.status = 1;
        res.x = detDx / detD;
        res.y = detDy / detD;
    } else if (detDx == 0 && detDy == 0)
        res.status = 0;
    else
        res.status = -1;
    return res;
}


string read_three_digits(int n) {
    const string a[10] = {"khong", "mot", "hai", "ba", "bon", "nam", "sau", "bay", "tam", "chin"};
    int tram = n / 100;
    int chuc = (n / 10) % 10;
    int donvi = n % 10;

    string result = "";
    if (tram > 0) {
        result += a[tram] + " tram ";
        if (chuc == 0 && donvi != 0)
            result += "le ";
    } else if (n > 0)
        result += "khong tram ";

    if (chuc > 1) {
        result += a[chuc] + " muoi ";
        if (donvi != 0)
            result += a[donvi] + " ";
    } else if (chuc == 1) {
        result += "muoi ";
        if (donvi != 0)
            result += a[donvi] + " ";
    } else if (chuc == 0 && donvi != 0 && tram != 0)
        result += a[donvi] + " ";
    else if (tram == 0 && chuc == 0 && donvi != 0)
        result += a[donvi] + " ";

    return result;
}

string read_number(int n) {
    if (n == 0)
        return "khong";

    const string unit[3] = {"", "nghin", "trieu"};
    int parts[3] = {0, 0, 0}; 


    parts[0] = n % 1000;      
    parts[1] = (n / 1000) % 1000; 
    parts[2] = n / 1000000;    

    string result = "";

    if (parts[2] != 0) {
        result += read_three_digits(parts[2]) + "trieu ";
    }

    if (parts[1] != 0) {
        result += read_three_digits(parts[1]) + "nghin ";
    } else if (parts[2] != 0 && (parts[0] != 0))
        result += "khong tram ";

    if (parts[0] != 0) {
        result += read_three_digits(parts[0]);
    }

    return format_read_number(result);
}

float convert_temperature(float temperature, bool f2c) {
    float result = f2c
        ? 5.0f * (temperature - 32) / 9.0f
        : 9.0f * temperature / 5.0f + 32;
    return round_to_decimal(result, 2);
}

float linear_regression(Point data_points[], int num_data_points, float new_data,
                        float learning_rate, int num_iterations) {
    float m = 0.0f, b = 0.0f;
    for (int s = 0; s < num_iterations; s++) {
        float delta_m = 0, delta_b = 0;
        for (int i = 0; i < num_data_points; i++) {
            float x = data_points[i].x;
            float y = data_points[i].y;
            delta_m += x * (y - (m * x + b));
            delta_b += (y - (m * x + b));
        }
        delta_m = (-2.0f / num_data_points) * delta_m;
        delta_b = (-2.0f / num_data_points) * delta_b;
        m -= learning_rate * delta_m;
        b -= learning_rate * delta_b;
    }
    return round_to_decimal(m * new_data + b, 3);
}
