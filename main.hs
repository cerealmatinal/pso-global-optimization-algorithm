module PSO where

import System.Random (getStdGen, randomRs)
import Data.List (minimumBy)
import Data.Function (on)

type Position = [Double]
type Velocity = [Double]
type Particle = (Position, Velocity)
type Swarm = [Particle]
type Bounds = [(Double, Double)]
type Objective = Position -> Double
type FitnessFunction = Position -> Double
type InertiaWeightFunction = Int -> Double

numDimensions :: Int
numDimensions = 2

numParticles :: Int
numParticles = 30

maxIterations :: Int
maxIterations = 100

bounds :: Bounds
bounds = replicate numDimensions (-5.0, 5.0)

c1 :: Double
c1 = 2.0

c2 :: Double
c2 = 2.0

inertiaWeight :: InertiaWeightFunction
inertiaWeight i = 0.9 - fromIntegral i * 0.5 / fromIntegral maxIterations

rosenbrock :: Objective
rosenbrock [x1, x2] = (1 - x1)^2 + 100 * (x2 - x1^2)^2

fitness :: FitnessFunction
fitness = negate . rosenbrock

initializeParticle :: StdGen -> Particle
initializeParticle gen = (position, velocity)
  where position = take numDimensions $ randomRs (l, u) gen
        velocity = take numDimensions $ repeat 0.0
        (l, u) = bounds !! 0

initializeSwarm :: StdGen -> Swarm
initializeSwarm gen = map initializeParticle $ take numParticles $ iterate snd $ split gen

updateVelocity :: Particle -> Position -> Velocity -> Velocity
updateVelocity (pos, vel) bestPos globalBest = newVelocity
  where newVelocity = zipWith (+) (zipWith (*) vel inertiaTerm) cognitiveTerm
        inertiaTerm = map (* inertia) vel
        cognitiveTerm = zipWith (*) (zipWith (-) bestPos pos) (repeat c1)
                     ++ zipWith (*) (zipWith (-) globalBest pos) (repeat c2)
        inertia = inertiaWeight iteration
        iteration = 0

updatePosition :: Particle -> Velocity -> Position
updatePosition (pos, vel) newVel = newPosition
  where newPosition = zipWith (+) pos newVel

updateParticle :: Particle -> Swarm -> Particle
updateParticle particle swarm = (newPosition, newVelocity)
  where bestParticle = bestPosition $ map fst swarm
        newVelocity = updateVelocity particle (fst particle) (fst bestParticle)
        newPosition = updatePosition particle newVelocity

updateSwarm :: Swarm -> Swarm
updateSwarm swarm = map (`updateParticle` swarm) swarm

bestPosition :: [Position] -> Position
bestPosition positions = positions !! bestIndex
  where values = map fitness positions
        bestIndex = fst $ minimumBy (compare `on` snd) $ zip [0..] values

bestParticle :: Swarm -> Particle
bestParticle swarm = fst $ minimumBy compareParticles swarm
  where
    compareParticles (_,p1) (_,p2) = compare (objective p1) (objective p2)

bestPosition :: Swarm -> Position
bestPosition swarm = position $ bestParticle swarm

updateVelocity :: RandomGen g => Parameters -> Particle -> Swarm -> g -> Velocity
updateVelocity params p swarm gen = newVelocity
  where
    c1 = c1Param params
    c2 = c2Param params
    r1 = fst $ randomR (0,1) gen
    r2 = fst $ randomR (0,1) gen
    oldVelocity = velocity p
    bestP = position $ bestParticle swarm
    bestS = bestPosition swarm
    newVelocity = zipWith (+) (zipWith (+) (map (c1 * r1 *) (zipWith (-) bestP (position p))) (map (c2 * r2 *) (zipWith (-) bestS (position p)))) (map (* w) oldVelocity)

updatePosition :: Particle -> Velocity -> Position
updatePosition p v = zipWith (+) (position p) v

updateParticle :: RandomGen g => Parameters -> Particle -> Swarm -> g -> Particle
updateParticle params p swarm gen = p { position = newPosition, velocity = newVelocity }
  where
    newPosition = updatePosition p newVelocity
    newVelocity = updateVelocity params p swarm gen

updateSwarm :: RandomGen g => Parameters -> Swarm -> g -> Swarm
updateSwarm params swarm gen = map (updateParticle params' swarm gen') swarm
  where
    params' = updateParameters params
    gen' = snd $ split gen

runPSO :: IO ()
runPSO = do
  gen <- getStdGen
  let swarm = initializeSwarm gen
      loop i bestSwarm
        | i >= maxIterations = return ()
        | otherwise = do
          let bestParticle = bestPosition swarm
              bestValue = objective bestParticle
              bestSwarm' = if bestValue < objective bestSwarm then bestParticle else bestSwarm
          putStrLn $ "Iteration " ++ show i ++ ": Best value = " ++ show bestValue
          loop (i + 1) bestSwarm'
            where
              swarm' = updateSwarm params swarm gen
              params = defaultParameters
  loop 0 (head swarm)