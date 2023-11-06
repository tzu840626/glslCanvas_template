#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;

float glow(float d, float str, float thickness){
    return thickness / pow(d, str);
}

vec2 hash2( vec2 x )            //亂數範圍 [-1,1]
{
    const vec2 k = vec2(-0.040,-0.680);
    x = x*k + k.yx;
    return -0.9 + 2.0*fract( 16.0 * k*fract( x.x*x.y*(x.x+x.y)) );
}
float gnoise( in vec2 p )       //亂數範圍 [-1,1]
{
    vec2 i = floor( p );
    vec2 f = fract( p );
    
    vec2 u = f*f*(3.0-2.0*f);

    return mix( mix( dot( hash2( i + vec2(0.0,0.0) ), f - vec2(0.0,0.0) ), 
                            dot( hash2( i + vec2(1.0,0.0) ), f - vec2(1.0,0.0) ), u.x),
                         mix( dot( hash2( i + vec2(0.0,1.0) ), f - vec2(0.0,1.0) ), 
                            dot( hash2( i + vec2(1.0,1.0) ), f - vec2(1.0,1.0) ), u.x), u.y);
}
#define Use_Perlin
//#define Use_Value
float noise( in vec2 p )        //亂數範圍 [-1,1]
{
#ifdef Use_Perlin    
return gnoise(p);   //gradient noise
#elif defined Use_Value
return vnoise(p);       //value noise
#endif    
return 0.0;
}
float fbm(in vec2 uv)       //亂數範圍 [-1,1]
{
    float f;                                                //fbm - fractal noise (4 octaves)
    mat2 m = mat2( 1.6,  1.2, -1.2,  1.6 );
    f   = 0.5000*noise( uv ); uv = m*uv;          
    f += 0.2500*noise( uv ); uv = m*uv;
    f += 0.1250*noise( uv ); uv = m*uv;
    f += 0.0625*noise( uv ); uv = m*uv;
    return f;
}

float sdStar5(in vec2 p, in float r, in float rf) //星星圖案
{
    const vec2 k1 = vec2(0.809016994375, -0.587785252292);
    const vec2 k2 = vec2(-k1.x,k1.y);
    p.x = abs(p.x);
    p -= 2.0*max(dot(k1,p),0.0)*k1;
    p -= 2.0*max(dot(k2,p),0.0)*k2;
    p.x = abs(p.x);
    p.y -= r;
    vec2 ba = rf*vec2(-k1.y,k1.x) - vec2(0,1);
    float h = clamp( dot(p,ba)/dot(ba,ba), 0.0, r );
    return length(p-ba*h) * sign(p.y*ba.x-p.x*ba.y);
}


///繪製愛心heart code
float M_SQRT_2=1.41421356237;
float heart(vec2 P, float size)
{
    float x = M_SQRT_2/2.0 * (P.x - P.y);
    float y = M_SQRT_2/2.0 * (P.x + P.y);
    float r1 = max(abs(x),abs(y))-size/3.5;
    float r2 = length(P - M_SQRT_2/2.0*vec2(+1.0,-1.0)*size/3.5)
    - size/3.5;
    float r3 = length(P - M_SQRT_2/2.0*vec2(-1.0,-1.0)*size/3.5)
    - size/3.5;
	return min(min(r1,r2),r3);
}

 
float sdCircleWave( in vec2 p, in float tb, in float ra )  //圖形 Circle Wave
{
    tb = 3.1415927*5.0/6.0*max(tb,0.0001);
    vec2 co = ra*vec2(sin(tb),cos(tb));
    p.x = abs(mod(p.x,co.x*4.0)-co.x*2.0);
    vec2  p1 = p;
    vec2  p2 = vec2(abs(p.x-2.0*co.x),-p.y+2.0*co.y);
    float d1 = ((co.y*p1.x>co.x*p1.y) ? length(p1-co) : abs(length(p1)-ra));
    float d2 = ((co.y*p2.x>co.x*p2.y) ? length(p2-co) : abs(length(p2)-ra));
    return min(d1, d2); 
}


float sdBlobbyCross( in vec2 pos, float he ) //圖形 lobbyCross
{
    pos = abs(pos);
    pos = vec2(abs(pos.x-pos.y),1.0-pos.x-pos.y)/sqrt(2.0);

    float p = (he-pos.y-0.25/he)/(6.0*he);
    float q = pos.x/(he*he*1.0);
    float h = q*q - p*p*p;
    
    float x;
    if( h>0.0 ) { float r = sqrt(h); x = pow(q+r,1.0/3.0)-pow(abs(q-r),1.0/3.0)*sign(r-q); }
    else        { float r = sqrt(p); x = 2.0*r*cos(acos(q/(p*r))/3.0); }
    x = min(x,sqrt(2.0)/2.0);
    
    vec2 z = vec2(x,he*(1.0-2.0*x*x)) - pos;
    return length(z) * sign(z.y);
    
    
}

float sdCross( in vec2 p, in vec2 b, float r ) 
{
    p = abs(p); p = (p.y>p.x) ? p.yx : p.xy;
    vec2  q = p - b;
    float k = max(q.y,q.x);
    vec2  w = (k>0.0) ? q : vec2(b.y-p.x,-k);
    return sign(k)*length(max(w,0.0)) + r;
}


//////////////////////////////////////////////////////


void main() {
    vec2 uv = gl_FragCoord.xy/u_resolution.xy;
    uv.x *= u_resolution.x/u_resolution.y;
    uv= uv*2.0-1.0;
    
    vec2 p = (2.0*uv-u_resolution.xy)/u_resolution.y;
    vec2 m = (2.0*u_mouse.xy-u_resolution.xy)/u_resolution.y;

    
    //陰晴圓缺
    float pi=3.14159;
    float theta=1.0*pi*u_time/9.0;
    vec2 point=vec2(sin(theta), cos(theta)); //隨著時間去做轉動
    float dir= dot(point, (uv))+0.5;
    
    //亂數作用雲霧
    float fog= fbm(0.4*uv+vec2(0.6*u_time, -0.01*u_time))*0.5+0.1;

    //定義圓環
    float dist = length(uv);
    float circle_dist = abs(dist-0.1);  //光環大小
    
    //定義迴圈
    float result;
    for(int index=0; index<7;++index) //迴圈線條的多寡
 {  
	//model
    vec2 uv_flip=vec2(-uv.x,-uv.y);	//宣告新的UV_flip 翻轉正向
    float weight=smoothstep(-0.1,0.6,-uv.y); 
    //float weight=step(uv.y,0.2); 	//(uv.y,0.2)顛倒 會交換上下noise的位置
    float freq=2.0+float(index)*1.0; //迴圈的線
 	float noise=gnoise(uv_flip*freq)*0.5*weight;
    float model_dist=abs(sdCircleWave(uv,1.0,0.4)+noise);//ads(sdStar5(uv,0.672,0.420))
        
        
        
 
    //動態呼吸
        
    // float breathing=sin(2.0*u_time/5.0*pi)*0.5+0.536;                     //option1
    float tb = 0.5-0.4*cos(u_time/0.9);
    // distance
    float d = sdCircleWave(p,tb, 0.5);

    float breathing=(exp(sin(u_time/5.0*pi)) - 0.36787944)*-0.2;         //option2 錯誤
     //float breathing=(exp(sin(u_time/2.0*pi)) - 0.36787944)*0.42545906412;                //option2 正確
    float strength =(0.3*tb+0.05);          //[0.2~0.3]         //光暈強度加上動態時間營造呼吸感
    float thickness=(0.2*breathing+0.1);          //[0.1~0.2]         //光環厚度 營造呼吸感
    // float glow_circle = glow(circle_dist, strength, thickness);
    // float glow_circle = glow(model_dist, strength, thickness);
    float glow_circle = glow(model_dist, strength, thickness);
    result+=glow_circle;
    // animation
    
    


}       
    gl_FragColor = vec4((vec3(result)+fog)*dir*vec3(0.658,0.8,1.000),2.0);
}