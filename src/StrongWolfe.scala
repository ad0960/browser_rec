package breeze.optimize

import breeze.util.SerializableLogging

/** 
 * Line Search方法
 */



abstract class CubicLineSearch extends SerializableLogging with MinimizingLineSearch {
  import scala.math._

  case class Bracket(
    t: Double, // 1d line search parameter
    dd: Double, // Directional Derivative at t
    fval: Double // Function value at t
    )

  /*
   * Invoke line search, returning stepsize
   */
  def minimize(f: DiffFunction[Double], init: Double = 1.0): Double
  /*
   * Cubic interpolation to find the minimum inside the bracket l and r.
   * Uses the fval and gradient at the left and right side, which gives
   * the four bits of information required to interpolate a cubic.
   * This is additionally "safe-guarded" whereby steps too close to
   * either side of the interval will not be selected.
   */
  def interp(l: Bracket, r: Bracket) = {
    // See N&W p57 actual for an explanation of the math
    val d1 = l.dd + r.dd - 3 * (l.fval - r.fval) / (l.t - r.t)
    val d2 = sqrt(d1 * d1 - l.dd * r.dd)
    val multipler = r.t - l.t
    val t = r.t - multipler * (r.dd + d2 - d1) / (r.dd - l.dd + 2 * d2)

    // If t is too close to either end bracket, move it closer to the middle

    val lbound = l.t + 0.1 * (r.t - l.t)
    val ubound = l.t + 0.9 * (r.t - l.t)
    t match {
      case _ if t < lbound =>
        logger.debug("Cubic " + t + " below LHS limit: " + lbound)
        lbound
      case _ if t > ubound =>
        logger.debug("Cubic " + t + " above RHS limit: " + ubound)
        ubound
      case _ => t
    }
  }
}

/*
 * StrongWolfe LineSearch方法
 */
class StrongWolfeLineSearch(maxZoomIter: Int, maxLineSearchIter: Int) extends CubicLineSearch {
  import scala.math._

  val c1 = 1e-4
  val c2 = 0.9

  def minimize(f: DiffFunction[Double], init: Double = 1.0):Double = {

    def phi(t: Double): Bracket = {
      val (pval, pdd) = f.calculate(t)
      Bracket(t = t, dd = pdd, fval = pval)
    }

    var t = init // Search's current multiple of pk
    var low = phi(0.0)
    val fval = low.fval
    val dd = low.dd

    if (dd > 0) {
      throw new FirstOrderException("Line search invoked with non-descent direction: " + dd)
    }

    def zoom(linit: Bracket, rinit: Bracket): Double = {

      var low = linit
      var hi = rinit

      for (i <- 0 until maxZoomIter) {
        // Interp assumes left less than right in t value, so flip if needed
        val t = if (low.t > hi.t) interp(hi, low) else interp(low, hi)

        // Evaluate objective at t, and build bracket
        val c = phi(t)
        //logger.debug("ZOOM:\n c: " + c + " \n l: " + low + " \nr: " + hi)
        logger.info("Line search t: " + t + " fval: " + c.fval +
          " rhs: " + (fval + c1 * c.t * dd) + " cdd: " + c.dd)


        if (c.fval > fval + c1 * c.t * dd || c.fval >= low.fval) {
          // "Sufficient decrease" condition not satisfied by c. Shrink interval at right
          hi = c
          logger.debug("hi=c")
        } else {

          // Zoom exit condition is the "curvature" condition
          // Essentially that the directional derivative is large enough
          if (abs(c.dd) <= c2 * abs(dd)) {
            return c.t
          }

          // If the signs don't coincide, flip left to right before updating l to c
          if (c.dd * (hi.t - low.t) >= 0) {
            logger.debug("flipping")
            hi = low
          }

          logger.debug("low=c")
          // If curvature condition not satisfied, move the left hand side of the
          // interval further away from t=0.
          low = c
        }
      }

      throw new FirstOrderException(s"Line search zoom failed")
    }

    ///////////////////////////////////////////////////////////////////

    for (i <- 0 until maxLineSearchIter) {
      val c = phi(t)

      // If phi has a bounded domain, inf or nan usually indicates we took
      // too large a step.
      if (java.lang.Double.isInfinite(c.fval) || java.lang.Double.isNaN(c.fval)) {
        t /= 2.0
        logger.error("Encountered bad values in function evaluation. Decreasing step size to " + t)
      } else {

        if ((c.fval > fval + c1 * t * dd) ||
          (c.fval >= low.fval && i > 0)) {
          logger.debug("Line search t: " + t + " fval: " + c.fval + " cdd: " + c.dd)
          return zoom(low, c)
        }

        if (abs(c.dd) <= c2 * abs(dd)) {
          return c.t
        }

        if (c.dd >= 0) {
          logger.debug("Line search t: " + t + " fval: " + c.fval +
            " rhs: " + (fval + c1 * t * dd) + " cdd: " + c.dd)
          return zoom(c, low)
        }

        low = c
        t *= 1.5
        logger.debug("Sufficent Decrease condition but not curvature condition satisfied. Increased t to: " + t)
      }
    }

    throw new FirstOrderException("Line search failed")
  }

}
