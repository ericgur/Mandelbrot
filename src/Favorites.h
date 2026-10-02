/**
 * @file Favorites.h
 * @brief Saved fractal locations ("favorites") and their storage in the application settings.
 *
 * A favorite records where the view is - center and zoom - together with the set type and
 * Julia constant that give those coordinates their meaning. Favorites are kept in QSettings,
 * so QCoreApplication's organization and application names must be set before they are
 * loaded or saved.
 */

#pragma once

#include <complex>
#include <optional>
#include <QRegularExpression>
#include <QString>
#include <QVector>
#include "QMandelbrotWidget.h"

/**
 * @struct Favorite
 * @brief A saved location in the fractal: view center, zoom, set type and Julia constant.
 *
 * The center keeps full fp128 precision. Beyond a zoom of about 2^44 a double can no longer
 * tell neighboring pixels apart, so a favorite rounded to double would land somewhere else.
 *
 * The Julia constant only matters for Julia favorites. It is stored for those alone, and
 * loading a Mandelbrot favorite leaves the current constant alone.
 */
struct Favorite {
    QString description;                                                           ///< Name shown in the Favorites menu.
    fp128_t centerX {};                                                            ///< Real part of the view center.
    fp128_t centerY {};                                                            ///< Imaginary part of the view center.
    int32_t log2Zoom = 0;                                                          ///< Zoom as a power of 2, in [logMinZoom, logMaxZoom].
    QMandelbrotWidget::set_type_t setType = QMandelbrotWidget::stMandelbrot;       ///< Fractal the location belongs to.
    std::complex<double> juliaConstant = QMandelbrotWidget::defaultJuliaConstant;  ///< Julia constant C; unused by Mandelbrot favorites.

    /**
     * @brief Get the name to show for this favorite.
     * @return The description, or a placeholder when it is blank.
     */
    [[nodiscard]] QString displayName() const;
};

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
 * @brief Parse a center coordinate typed or stored as decimal text.
 * @param text Text that must fully match CoordinateRegularExpression().
 * @return The value, accurate to about 36 decimal digits, or no value if the text doesn't match.
 */
[[nodiscard]] std::optional<fp128_t> ParseCoordinate(const QString& text);

/**
 * @brief Format a center coordinate as decimal text that ParseCoordinate() reads back.
 *
 * The text is the shortest decimal that ParseCoordinate() reads back as exactly the same
 * value, so "-0.2" stays "-0.2" although 0.2 has no exact binary form. When no shorter text
 * reads back exactly, every meaningful digit is written, which reads back to within one unit
 * in the last place.
 *
 * @param value Coordinate to format.
 * @return Decimal text, for example "-0.7436438870371587047521915034129".
 */
[[nodiscard]] QString CoordinateToString(const fp128_t& value);

/**
 * @brief Read the favorites from the application settings.
 *
 * Entries that can't be placed, such as ones with a coordinate that doesn't parse, are
 * skipped. That only happens when the settings were edited outside the application.
 * Because they are skipped, the next SaveFavorites() call drops them for good.
 *
 * @return The favorites in menu order; empty if none were saved.
 */
[[nodiscard]] QVector<Favorite> LoadFavorites();

/**
 * @brief Replace the favorites in the application settings.
 * @param favorites The complete list, in menu order.
 * @return True if the settings were written.
 */
bool SaveFavorites(const QVector<Favorite>& favorites);
