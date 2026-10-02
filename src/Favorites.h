/**
 * @file Favorites.h
 * @brief Saved fractal locations ("favorites") and their storage in the application settings.
 *
 * A favorite records where the view is - center and zoom - together with the set type and
 * Julia constant that give those coordinates their meaning, and the iteration limit the view
 * was rendered with. Favorites are kept in QSettings,
 * so QCoreApplication's organization and application names must be set before they are
 * loaded or saved.
 */

#pragma once

#include <QString>
#include <QVector>
#include "Fp128Text.h"
#include "QMandelbrotWidget.h"

/**
 * @struct Favorite
 * @brief A saved location in the fractal: view center, zoom, iteration limit, set type and Julia constant.
 *
 * The center keeps full fp128 precision. Beyond a zoom of about 2^44 a double can no longer
 * tell neighboring pixels apart, so a favorite rounded to double would land somewhere else.
 *
 * The iteration limit is part of what the view looked like: it decides how much of the
 * boundary resolves and how the palette spreads over it. Auto is stored as Auto rather than
 * as the limit it picked, which at the favorite's zoom is the same limit, and which keeps
 * scaling if the user zooms on from there.
 *
 * The Julia constant only matters for Julia favorites. It is stored for those alone, and
 * loading a Mandelbrot favorite leaves the current constant alone. Like the center, it keeps
 * full fp128 precision.
 */
struct Favorite {
    QString description;                                         ///< Name shown in the Favorites menu.
    fp128_t centerX {};                                          ///< Real part of the view center.
    fp128_t centerY {};                                          ///< Imaginary part of the view center.
    int32_t log2Zoom = 0;                                        ///< Zoom as a power of 2, in [logMinZoom, logMaxZoom].
    int64_t maxIterations = QMandelbrotWidget::auto_iterations;  ///< Iteration limit in [min_fixed_iterations, max_iterations], or auto_iterations.
    QMandelbrotWidget::set_type_t setType = QMandelbrotWidget::stMandelbrot;  ///< Fractal the location belongs to.
    Complex128 juliaConstant = QMandelbrotWidget::defaultJuliaConstant();     ///< Julia constant C; unused by Mandelbrot favorites.

    /**
     * @brief Get the name to show for this favorite.
     * @return The description, or a placeholder when it is blank.
     */
    [[nodiscard]] QString displayName() const;
};

/**
 * @brief Read the favorites from the application settings.
 *
 * When the settings hold no favorites list at all, as on a first run, this returns ten
 * famous Mandelbrot and Julia locations instead, none of them at zoom 1x. They reach the
 * settings with the first save. A list the user emptied is stored as empty and stays empty.
 *
 * Entries that can't be placed, such as ones with a coordinate that doesn't parse, are
 * skipped. That only happens when the settings were edited outside the application.
 * Because they are skipped, the next SaveFavorites() call drops them for good.
 *
 * @return The favorites in menu order; the defaults if no list was ever saved.
 */
[[nodiscard]] QVector<Favorite> LoadFavorites();

/**
 * @brief Replace the favorites in the application settings.
 * @param favorites The complete list, in menu order.
 * @return True if the settings were written.
 */
bool SaveFavorites(const QVector<Favorite>& favorites);
