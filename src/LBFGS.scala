package breeze.optimize

import breeze.linalg._
import breeze.linalg.operators.OpMulMatrix
import breeze.math.MutableInnerProductModule
import breeze.optimize.linear.PowerMethod
import breeze.util.SerializableLogging

/**
 * LBFGS 算法
 */
class LBFGS[T](maxIter: Int = -1, m: Int=10, tolerance: Double=1E-9)
              (implicit space: MutableInnerProductModule[T, Double]) extends FirstOrderMinimizer[T, DiffFunction[T]](maxIter, tolerance, tolerance) with SerializableLogging {

  import space._
  require(m > 0)

  type History = LBFGS.ApproximateInverseHessian[T]


  override protected def adjustFunction(f: DiffFunction[T]): DiffFunction[T] = f.cached

  protected def takeStep(state: State, dir: T, stepSize: Double) = state.x + dir * stepSize
  protected def initialHistory(f: DiffFunction[T], x: T):History = new LBFGS.ApproximateInverseHessian(m)
  protected def chooseDescentDirection(state: State, fn: DiffFunction[T]):T = {
    state.history * state.grad
  }

  protected def updateHistory(newX: T, newGrad: T, newVal: Double,  f: DiffFunction[T], oldState: State): History = {
    oldState.history.updated(newX - oldState.x, newGrad :- oldState.grad)
  }

  /**
   * @param state the current state
   * @param f The objective
   * @param dir The step direction
   * @return stepSize
   */
  protected def determineStepSize(state: State, f: DiffFunction[T], dir: T) = {
    val x = state.x
    val grad = state.grad

    val ff = LineSearch.functionFromSearchDirection(f, x, dir)
    val search = new StrongWolfeLineSearch(maxZoomIter = 10, maxLineSearchIter = 10) // TODO: Need good default values here.
    val alpha = search.minimize(ff, if(state.iter == 0.0) 1.0/norm(dir) else 1.0)

    if(alpha * norm(grad) < 1E-10)
      throw new StepSizeUnderflow
    alpha
  }
}

object LBFGS {
  case class ApproximateInverseHessian[T](m: Int,
                                          private[LBFGS] val memStep: IndexedSeq[T] = IndexedSeq.empty,
                                          private[LBFGS] val memGradDelta: IndexedSeq[T] = IndexedSeq.empty)
                                         (implicit space: MutableInnerProductModule[T, Double]) extends NumericOps[ApproximateInverseHessian[T]] {

    import space._

    def repr: ApproximateInverseHessian[T] = this

    def updated(step: T, gradDelta: T) = {
      val memStep = (step +: this.memStep) take m
      val memGradDelta = (gradDelta +: this.memGradDelta) take m

      new ApproximateInverseHessian(m, memStep,memGradDelta)
    }


    def historyLength = memStep.length

    def *(grad: T) = {
     val diag = if(historyLength > 0) {
       val prevStep = memStep.head
       val prevGradStep = memGradDelta.head
       val sy = prevStep dot prevGradStep
       val yy = prevGradStep dot prevGradStep
       if(sy < 0 || sy.isNaN) throw new NaNHistory
       sy/yy
     } else {
       1.0
     }

     val dir = space.copy(grad)
     val as = new Array[Double](m)
     val rho = new Array[Double](m)

     for(i <- 0 until historyLength) {
       rho(i) = (memStep(i) dot memGradDelta(i))
       as(i) = (memStep(i) dot dir)/rho(i)
       if(as(i).isNaN) {
         throw new NaNHistory
       }
       axpy(-as(i), memGradDelta(i), dir)
     }

     dir *= diag

     for(i <- (historyLength - 1) to 0 by (-1)) {
       val beta = (memGradDelta(i) dot dir)/rho(i)
       axpy(as(i) - beta, memStep(i), dir)
     }

     dir *= -1.0
     dir
    }
  }

  implicit def multiplyInverseHessian[T](implicit vspace: MutableInnerProductModule[T, Double]):OpMulMatrix.Impl2[ApproximateInverseHessian[T], T, T] = {
    new OpMulMatrix.Impl2[ApproximateInverseHessian[T], T, T] {
      def apply(a: ApproximateInverseHessian[T], b: T): T = a * b
    }
  }
}

