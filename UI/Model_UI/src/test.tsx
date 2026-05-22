import { useMemo } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer } from 'recharts';

interface args {
    x: Array<number>,
    y: Array<number>,
    name: string
}

const MatplotLibStyleChart = ({ x, y , name}: args) => {
  // 1. Optimization: Downsample and transform data inside useMemo
  // This prevents 32,000 objects being recreated and re-diffed on every render.
  const data = useMemo(() => {
    const totalPoints = x.length;
    const maxVisiblePoints = 2000; // Optimal for performance vs visual detail
    const step = Math.ceil(totalPoints / maxVisiblePoints);

    const result = [];
    for (let i = 0; i < totalPoints; i += step) {
      result.push({
        name: x[i],
        value: y[i],
      });
    }
    return result;
  }, [x, y]);

  return (
    <div style={{ width: '100%', height: 400, backgroundColor: '#f9f9f9', padding: '20px' }}>
      <h3>{name}</h3>
      <ResponsiveContainer width="100%" height="100%">
        <LineChart data={data}>
          <CartesianGrid strokeDasharray="3 3" />
          <XAxis 
            dataKey="name" 
            label={{ value: 'X Axis', position: 'insideBottom', offset: -5 }}
            minTickGap={30} // Prevents X-Axis labels from overlapping and lagging
          />
          <YAxis label={{ value: 'Y Axis', angle: -90, position: 'insideLeft' }} />
          
          {/* 2. Optimization: Disable tooltip animation to prevent lag on mouse move */}
          <Tooltip isAnimationActive={false} />
          
          <Line 
            type="linear"      // "linear" is much faster to calculate than "monotone"
            dataKey="value" 
            stroke="#8884d8" 
            strokeWidth={2}
            dot={false}        // CRITICAL: Rendering 32,000 <circle> nodes will crash the browser
            activeDot={{ r: 6 }} 
            isAnimationActive={false} // CRITICAL: Prevents heavy CPU usage on load
          />
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
};

export default MatplotLibStyleChart;