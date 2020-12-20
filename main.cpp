#include <algorithm>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>

struct Point {
    double x;
    double y;
};

struct Coefficients {
    double c;
    double d;
};

double linearFunction(const double _x) {
    return -0.5*_x + 0;
}

double error(const std::vector<Point>& right, const double c, const double d) {
    double sum = 0.;
    for (auto point : right) {
        sum += pow(point.y - (c * point.x + d), 2);
    }
    return sum;
}

double randomInRange(const double lower, const double upper) {
    return lower + rand() * 1./RAND_MAX * (upper - lower);
}

std::vector<Point> divideInterval(const double lower, const double upper,
                                  const size_t pointsNumber, const double fault) {
    std::vector<Point> points(pointsNumber);
    const double step = (upper - lower) / static_cast<double>(pointsNumber - 1);
    for (size_t i = 0; i < pointsNumber; ++i) {
        points[i].x = lower + i * step;
        points[i].y = linearFunction(points[i].x) + randomInRange(- fault / 2, fault / 2);
    }
    return points;
}

void findBestDCoefficientForC(const std::vector<Point>& points, Coefficients& current) {
    const double D_MIN = -2;
    const double D_MAX = 2;
    const size_t iterations = 50;

    current.d = D_MIN;
    for (size_t i = 0; i < iterations; ++i) {
        double new_d = randomInRange(D_MIN, D_MAX);
        if (error(points, current.c, new_d) < error(points, current.c, current.d)) {
            current.d = new_d;
        }
    }
}

Coefficients findBestCoefficients(const std::vector<Point>& points) {
    const double C_MIN = -2.5;
    const double C_MAX = 1.5;
    const size_t iterations = 50;
    const double step = (C_MAX - C_MIN) / static_cast<double>(iterations - 1);

    std::vector<Coefficients> allVariants;
    for (size_t i = 0; i < iterations; ++i) {
        Coefficients newCoefficients{};
        newCoefficients.c = C_MIN + i * step;
        findBestDCoefficientForC(points, newCoefficients);
        allVariants.push_back(newCoefficients);
    }

    Coefficients bestCoefficients = allVariants[0];
    for (const Coefficients& item : allVariants) {
        if (error(points, item.c, item.d) < error(points, bestCoefficients.c, bestCoefficients.d)) {
            bestCoefficients = item;
        }
    }
    return bestCoefficients;
}

void printTable(const std::vector<Point>& points) {
    std::cout << std::string(27, '-') << std::endl;
    std::cout << "| " << std::left << std::setw(10) << 'x'
              << " | " << std::setw(10) << "f(x)" << " |\n";
    std::cout << std::string(27, '-') << std::endl;
    for (const auto& item : points) {
        std::cout << "| " << std::setw(10) << item.x
                  << " | " << std::setw(10) << item.y << " |\n";
    }
    std::cout << std::string(27, '-') << std::endl;
}

int main() {
    const double MINIMUM = -2.;
    const double MAXIMUM = 2.;
    const size_t POINTS = 16;
    const double ERROR_LIMIT = 2.;

    // Part 1. Function without noise
    srand(time(nullptr));
    auto points_1 = divideInterval(MINIMUM, MAXIMUM, POINTS, 0.);
    auto coefficients_1 = findBestCoefficients(points_1);
    std::cout << "Part 1. Without noise\n"
              << "Right function: y = -0.5x\nTable of values:\n";
    printTable(points_1);
    std::cout << "Found function: y = " << coefficients_1.c
              << "x " << (coefficients_1.d >= 0 ? "+ ": "- ")
              << (coefficients_1.d >= 0 ? coefficients_1.d : (coefficients_1.d * -1.)) << std::endl;

    // Part 2. Function with noise
    auto points_2 = divideInterval(MINIMUM, MAXIMUM, POINTS, ERROR_LIMIT);
    auto coefficients_2 = findBestCoefficients(points_2);
    std::cout << "Part 2. With noise\n"
              << "Right function: y = -0.5x + random(A)\nTable of values:\n";
    printTable(points_2);
    std::cout << "Found function: y = " << coefficients_2.c
              << "x " << (coefficients_2.d >= 0 ? "+ ": "- ")
              << (coefficients_2.d >= 0 ? coefficients_2.d : (coefficients_2.d * -1.)) << std::endl;

    return 0;
}
