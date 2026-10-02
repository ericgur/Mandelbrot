/**
 * @file Fp128Text.cpp
 * @brief Implementation of the decimal text conversions for fp128 values.
 */

#include "pch.h"
#include <string>
#include "Fp128Text.h"

using namespace Qt::StringLiterals;

namespace
{

/**
 * @brief Round a decimal number's text to fewer fraction digits, half up.
 * @param text Decimal number with a fraction, as fp128 formats it, for example "-0.25".
 * @param end Index of the first digit to drop; it must lie past the decimal point.
 * @return The text up to @p end, rounded, for example "-0.3" for "-0.25" and an @p end of 4.
 */
[[nodiscard]] std::string RoundDecimal(const std::string& text, size_t end)
{
    std::string result = text.substr(0, end);
    if (text[end] < '5') {
        return result;
    }

    // carry through the kept digits, stepping over the decimal point
    for (size_t i = result.size(); i-- > 0;) {
        char& c = result[i];
        if (c == '-') {
            break;
        }
        if (c == '.') {
            continue;
        }
        if (c != '9') {
            ++c;
            return result;
        }
        c = '0';
    }
    // every kept digit was a 9, so the carry adds a digit: "9.96" -> "10.0"
    result.insert(result.starts_with('-') ? 1 : 0, 1, '1');

    return result;
}

/**
 * @brief Parse decimal text that must match a pattern in full.
 * @param anchored Pattern anchored at both ends.
 * @param text Text to parse.
 * @return The value, or no value if the text doesn't match.
 */
[[nodiscard]] std::optional<fp128_t> ParseMatching(const QRegularExpression& anchored, const QString& text)
{
    if (!anchored.match(text).hasMatch()) {
        return std::nullopt;
    }

    return fp128_t(text.toStdString());
}

}  // namespace

const QRegularExpression& CoordinateRegularExpression()
{
    static const QRegularExpression pattern(uR"([+-]?(\d{1,2}(\.\d*)?|\.\d+))"_s);
    return pattern;
}

const QRegularExpression& JuliaComponentRegularExpression()
{
    // 0 and 1 take any fraction, 2 only zeros
    static const QRegularExpression pattern(uR"([+-]?([01](\.\d*)?|2(\.0*)?|\.\d+))"_s);
    return pattern;
}

std::optional<fp128_t> ParseCoordinate(const QString& text)
{
    static const QRegularExpression anchored(QRegularExpression::anchoredPattern(CoordinateRegularExpression().pattern()));
    return ParseMatching(anchored, text);
}

std::optional<fp128_t> ParseJuliaComponent(const QString& text)
{
    static const QRegularExpression anchored(QRegularExpression::anchoredPattern(JuliaComponentRegularExpression().pattern()));
    return ParseMatching(anchored, text);
}

QString ToDecimalText(const fp128_t& value)
{
    // fp128 prints every meaningful digit, so a value typed as "-0.2" would come back as
    // "-0.200000000000000000000000000000000001", 0.2 having no exact binary form. The fewest
    // digits that parse back to the same value keep such numbers as they were typed, and
    // checking each candidate with the parser also makes the round trip exact where possible.
    const std::string full = static_cast<std::string>(value);
    const size_t point = full.find('.');
    if (point == std::string::npos) {
        return QString::fromStdString(full);
    }

    for (size_t end = point + 2; end < full.size(); ++end) {
        const std::string candidate = RoundDecimal(full, end);
        if (fp128_t(candidate) == value) {
            return QString::fromStdString(candidate);
        }
    }

    return QString::fromStdString(full);
}
