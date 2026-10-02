/**
 * @file Fp128Text.h
 * @brief Decimal text for the fp128 values the UI shows, takes as input and stores.
 *
 * Text is the only form that carries an fp128 value's full precision through a line edit or
 * the settings store; a double keeps about 17 significant digits of the 36 an fp128 value
 * holds. Center coordinates and Julia constant parts both go through these functions, each
 * with its own range.
 */

#pragma once

#include <optional>
#include <QRegularExpression>
#include <QString>
#include "QMandelbrotWidget.h"

/**
 * @brief Get the regular expression a center coordinate's text must match.
 *
 * Accepts a plain decimal number with an optional sign and at most two integer digits, for
 * example "-0.7436438870371587047521915034129". Exponents are not accepted. The two-digit
 * limit keeps a coordinate, plus the view's half-width around it, well inside the
 * [-128, 128) range of fp128_t.
 *
 * The pattern is not anchored, which suits QRegularExpressionValidator: it anchors the
 * pattern itself and reports a partial match, such as a lone "-", as Intermediate.
 *
 * @return The coordinate pattern.
 */
[[nodiscard]] const QRegularExpression& CoordinateRegularExpression();

/**
 * @brief Get the regular expression the text of a Julia constant's real or imaginary part must match.
 *
 * Accepts a plain decimal number in [-2, 2], the range in which Julia sets are not just dust,
 * for example "-0.122561166876653619975245551820735654". Like CoordinateRegularExpression(),
 * it takes no exponents and is not anchored.
 *
 * @return The Julia constant part pattern.
 */
[[nodiscard]] const QRegularExpression& JuliaComponentRegularExpression();

/**
 * @brief Parse a center coordinate typed or stored as decimal text.
 * @param text Text that must fully match CoordinateRegularExpression().
 * @return The value, accurate to about 36 decimal digits, or no value if the text doesn't match.
 */
[[nodiscard]] std::optional<fp128_t> ParseCoordinate(const QString& text);

/**
 * @brief Parse a Julia constant's real or imaginary part typed or stored as decimal text.
 * @param text Text that must fully match JuliaComponentRegularExpression().
 * @return The value, accurate to about 36 decimal digits, or no value if the text doesn't match.
 */
[[nodiscard]] std::optional<fp128_t> ParseJuliaComponent(const QString& text);

/**
 * @brief Format an fp128 value as decimal text that the parse functions read back.
 *
 * The text is the shortest decimal that reads back as exactly the same value, so "-0.2" stays
 * "-0.2" although 0.2 has no exact binary form. When no shorter text reads back exactly,
 * every meaningful digit is written, which reads back to within one unit in the last place.
 *
 * @param value Value to format.
 * @return Decimal text, for example "-0.7436438870371587047521915034129".
 */
[[nodiscard]] QString ToDecimalText(const fp128_t& value);
